# Test Architect — Lessons (hot digest)

Read at the START of every dispatch. This is the **behavioral index**: one rule
per entry, imperative and standalone. War stories, measured numbers and
`file:line` detail live in **`lessons_archive.md`** (sections L1–L61) — open only
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
- **⛔⛔ RETIRING A DEFERRAL: never read the reason string — RUN the row. The fix
  the exclusion CITES is usually not the fix that heals it.** `[M]` 4 of 5
  imperative `pytest.xfail` sites all named issue #200 as their re-enabler; #200
  is still OPEN, yet all four PASS today (`rel` 3.6e-15 / 4.5e-13 / 1.2e-14 vs
  the closed form) — cured by an unrelated `restart`-sizing lineage (ERR-053, then
  #282). Checking "has the cited issue landed?" answers "still blocked" and is
  WRONG. Method: a `-p` plugin monkeypatching `pytest.xfail` to a no-op that
  ASSERTS its own installation (a `sessionfinish` neutralisation COUNT), then one
  run per row. Attribute the healing afterwards by grepping the production
  lineage — and if two candidate fixes exist, say the cure is in the lineage and
  that the decisive one is UNDISCRIMINATED rather than picking one. Corollaries:
  ⭐ a reason string can be false in a SECOND way nobody noticed — the budget it
  names may not be a live knob (`max_inner` → scipy `maxiter` = restart CYCLES
  with `restart == n_dof`, so `max_inner=2` returns the `max_inner=1000` answer to
  the last bit); ⛔ an UNCONDITIONAL stub whose body is ONLY the `pytest.xfail`
  call cannot be converted literally — a strict marker over an empty body XPASSes
  and reds, so the conversion must SUPPLY a body (a concept-level capability
  probe: grep the *words* the issue uses, not one guessed symbol) where exactly one
  statement can fail and it is the documented reason; ⛔ the healed row's doc
  claims are in the blast radius — a page calling it "the xfail'd cell" and naming
  the file goes present-tense-FALSE (`vv` anti-#21). → `L45`
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
- **⛔ A new guard wired AFTER an existing one: the earlier guard's inputs are a
  DISCRIMINATION row, not a negative row.** `[M]` slab GL on a cylinder is
  classified `on_edge_node=True` by the predicate (it does NOT raise) yet is
  intercepted by the older structure-less guard, so a row asserting the NEW
  message is a false red and one asserting the OLD message twins two committed
  gates. The honest row asserts, on that same input, **old fragment PRESENT +
  new fragment ABSENT** — which pins the WIRING ORDER (move the helper one line
  earlier ⟹ reds). Mint the two messages with disjoint fragments and assert the
  disjointness once. → `L43c`
- **⛔ A ∀-quantifier over a per-element predicate is UNGATE-ABLE when no factory
  produces a MIXED input.** `[M]` every refused quadrature family had ALL levels
  non-carrying and every admitted one ALL carrying ⟹ `all(...)` and `any(...)`
  agree on every constructible input (Mode-12 at the FIXTURE). The fix is
  architectural: split the pure predicate out so it takes the element tuple and
  **returns the offending POSITIONS** (`vv` anti-#14), making a synthetic mixed
  tuple a two-line unit test. Without the split the quantifier ships provably
  ungated and the implementer is credited for it. → `L43d`
- **⭐ Count the rows that REACH the assertion, not the rows that exist.** A
  guard-clause early return (`if not seed_levels: check_absence(); return`) made
  `[M]` **10 of 20** rows of a battery return after one check — the invariant the
  test is NAMED for had only ever run on the other 6, each at a degenerate
  one-element size. Anti-#20 row inflation wearing a guard clause; read the
  early-exit branch before believing a parametrize count. → `L43f`
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
- **⛔⛔ A `warnings.warn(stacklevel=N)` is a claim about EVERY call site's
  DEPTH, and NO message gate can see it — the message is a pure function of
  its arguments, the attribution is a second observable.** `[M]` #340 N4.7:
  `stacklevel=3` was already false at 2 of 8 sites (two calls sat in PRIVATE
  helpers one frame below the public entry), so the warning blamed
  `orpheus/sn/solver.py:3541` — the library's own dispatch line — at every
  budget on two fixtures; `grep -rn stacklevel tests/` = **1 hit**, unrelated.
  ⭐⭐ **And the OBVIOUS gate is Mode-12 blind to half the class**: "the
  attributed file is not under the package" reds for `stacklevel→2` (blames
  the entry) and stays **GREEN** for `stacklevel→4` (blames the caller's
  caller — still outside the package, still the wrong line). Ship TWO legs:
  the portable `not is_relative_to(pkg_root)` over every entry, PLUS
  `w.lineno == inspect.currentframe().f_lineno + 1` recorded immediately
  above a DIRECT call. Never gate by reading the literal `stacklevel == 3` —
  signature-tautological, since the hazard is the call DEPTH, not the value.
  Structural fix beats gate: hoist every emission into the public entry so
  the depth is uniform. → `L46b`, `L46c`
- **⛔ A guard WIDENING (top-level → whole-tree predicate) needs a
  CHILD-failing fixture per family, or the battery is blind to the widening
  it is credited with.** `[M]` re-installing the pre-widening guard
  partitions exactly: child-starved rows go silent, outer-starved rows still
  fire. And a family whose inner is DIRECT (diffusion's LU, `budget == 0`)
  structurally cannot supply one — say so rather than counting its row
  (anti-#20). → `L46d`
- **⛔⛔ A "functional X is BLIND to `ker A`" gate whose null vector comes from
  `svd(A)` is a TAUTOLOGY — `A·(a null vector of A) = 0` is a fact about the
  FACTORISATION, measurable with no solver in the room.** Both blindness legs
  of the promoted #344 diagnostic were of this shape; only its two anchor legs
  and its control had teeth, and the docstring credited the tautological ones.
  ⭐ Re-pose onto the PRODUCTION stopping path: run the driver TWICE from cold
  starts differing only inside `ker A` (`v` is a FIXED direction of
  `G = M⁻¹N`, so SI preserves it exactly). `[M]` iteration count **344 both**,
  residual `9.028098e-14` vs `9.022488e-14`, balance projection `2.795085`
  BIT-identically, bulk `8.4e-16` — traces **11.26 %** apart, difference equal
  to `v` to `2.3e-14`; kernel-free control `2.6e-15`. Test the MEASURAND the
  production code reports, never a re-derivation of it. → `L49a`
- **⛔⛔ A RATIO-valued pin is Mode-12 blind to a uniform scale — write the
  measurand's stabiliser down BEFORE the row, then pin one UN-normalised quantity
  from the same object.** `[M]` CS3: multiplying EVERY norm by `(1+1e-9)` left the
  whole 11-point ρ trajectory GREEN (ρ is a ratio; a common factor cancels exactly)
  and reddened only the separate `‖Δψ‖` pin — so a gate pinning "the ρ trajectory"
  alone cannot see any uniform mis-scaling, including "the relocation forgot a
  weight". Same family: normalised shapes annihilate global scaling, spectra
  annihilate similarity. ⭐ Companion: also RUN the mutation you expect to be GREEN
  (`Field.l2 → np.linalg.norm`, `[M]` 5 passed) — a declared blindness never
  executed is an unmeasured claim wearing a measurement's clothes. → `L58b`
- **⭐ Before minting a bit-identity instrument, grep for a WARNING class you can
  ESCALATE — a project with an audible-drift tripwire has already built the gate.**
  `[M]` `-W error::tests.sn.regression._regression_assert.DriftWarning` turns the
  #208/#333 stored-value gate into a **1-ULP** wall on 3 drivers in **1.60 s**
  (control: a plugin advancing the first element of every loaded baseline by 1 ULP
  → 3 failed), and the same flag shows **11 of 13** DD regression cases bit-exact
  with the 2 exceptions NAMED. Verify both that the `-W` string PARSES (`vv` Mode-8
  EIGHTH class) and that it bites. → `L58c`
- **⭐⭐ MEASURE a vv#19 CONTROL on the pre-carve tree BEFORE designing the pair —
  a control can sit in the error's STABILISER, and both readings are small
  residuals.** `vv` #12 gives the OPERATOR-side blindness (`[G, Aᵀ]=0`); the dual
  nobody had written is the **SPACE-side**: a rank-1 (point) axis makes EVERY
  weight a one-element array = a SCALAR, and a scalar `G` commutes with
  everything. `[M]` on the homogeneous `(ng,1)` carrier, `loss.H.apply([[1],[2]])`
  reads `[-0.08, 0.2]` for `None`, `w=1` AND **`w=2`** (the design's proposed
  control — a provable non-catcher) and only `[-0.38, 0.2]` for a per-GROUP
  `w=[2,5]`. ⟹ the control had to become a deliberately NON-PHYSICAL toy, said so
  in its docstring. ⚠ Only component 0 moves — assert the whole vector. And the
  corollary: **when a distinction is invisible to every VALUE functional the
  carrier admits, IDENTITY is the only instrument** (here: the derived space NAME
  must encode the weights), so the battery needs a MUST-STAY-GREEN column naming
  the value gates that provably cannot see it. → `L59a`
- **⭐ A TYPE-ANNOTATION widen has NO runtime witness — its gate is the type
  checker.** "Hand it the wider type, assert it constructs" is green BEFORE and
  after (Python does not enforce annotations): `vv` Mode-8's signature-tautological
  class wearing a type hint. Gate it with `tests/test_pyright_ratchet.py` and record
  the delta in the commit. Same pass, same family: deleting a leftover
  `basis_shape=(ng,1)` once the domain derives it is **value-identical**, so
  "both spellings gone" is a **grep obligation**, not a gate — say which done-when
  items are grep obligations, or they read as covered. → `L59d`
- **⛔ An attribute→property conversion kills every committed `hasattr(Class, …)`
  PREMISE — grep them BEFORE the carve, not after the red.** `[M]`
  `tests/test_docstring_xrefs.py:391` asserts `not hasattr(SNMesh, "mesh")` as the
  *premise* of the unannotated-instance-attribute row; a forwarding property on the
  base makes it True and the row reds on its premise line, not on its subject. A
  600-line census missed it because it swept `.mesh` READERS and a premise is not a
  read. The replacement is one `ast` scan away — `[M]` the hierarchy has exactly TWO
  bare public `self.X = …` attributes (`mesh`, `quad`) and `quad` satisfies all three
  assertions untouched. `grep -rn 'hasattr(' tests/` filtered to the converted names
  cost one command and returned one hit for ten names. → `L60c`
- **⛔ A guard placed after the CALLER's own attribute reads is unreachable — and the
  error TYPE is the placement gate.** `SNMesh.from_material_mesh` reads
  `material_mesh.axes`/`.mesh` in its own body, so a typed refusal inside `_init_core`
  never runs: the forwarding property's `AttributeError` fires a frame earlier. Put
  the guard at the TOP of the promoting classmethod, and note that
  `pytest.raises(ValueError, match=…)` then covers BOTH "it refuses" and "in the right
  place", since `AttributeError` is not a `ValueError`. ⭐ Companion: deleting the new
  guard still reds — **by raising** — so only the `match=` leg attributes it (L31).
  → `L60f`

- **⛔⛔ A charter's flagship NUMERICAL gate can be a THEOREM with no reachable
  falsifier — and the tell is in its own justification: a subordinate clause
  saying "…which today produces the same result for the WRONG reason".** `[M]`
  CS4a: `.H == apply_transpose` reads `0.000e+00` under the R2 defect
  (`space=None`) AND under the fix (`space=` the quotient), and `≤2.2e-16` on a
  meshed bulk of **56 000×** volume spread. Closed form: `A†=G⁻¹AᵀG`, all four
  leaves are spatially diagonal, every reachable bulk metric is
  `V_cell ⊗ counting` ⟹ `[G,Aᵀ]=0`; the only loading axis is a per-GROUP energy
  weight, which `EnergyAxis` **refuses at construction** — a construction
  invariant the SAME campaign shipped one phase earlier. Repair: gate the
  theorem's **PREMISE** (`apply_metric(x)` is `array_equal` to `x`; both axes
  `weights is None`) which IS red-capable, keep ONE corollary row labelled
  claim-kind THEOREM carrying the blindness table, and name the pre-existing
  vv#19 control as the only loaded partner. ⛔ Do NOT manufacture a wrong-metric
  control on the production mint — the type refuses it, so the control is
  unconstructible. **Grep a charter for "for the wrong reason" / "the same
  result today".** → `L61a`
- **⛔ A marker SPLIT (function-level `@xfail` → per-row `pytest.param(marks=)`)
  fails SILENTLY by losing `strict`, and both obvious gates are blind.** `[M]`
  `--collect-only` node-ids (98 lines) and `-rx` status+reason (16 lines) are
  BOTH unchanged when the new mark is spelled `pytest.mark.xfail(reason=…)`
  without `strict=True` — the row still reports `x`, the suite is green, and
  `pyproject.toml` has **no `xfail_strict`** so the default is non-strict. The
  catcher is 5 lines and permanent: `pytest.param(...).marks` is introspectable
  at import (`[('xfail', True)]` vs `[]`), so assert `strict is True` on every
  xfail mark of every row constant. ⭐ A *dropped* mark is by contrast a visible
  RED (the unmarked row fails on its own subject) — only the `strict` half is
  silent. → `L61g`

- **⭐⭐ When the phase's claim is "consumer X now reads from OWNER A instead of
  OWNER B", no value/identity gate states it — build a ROUTE gate: pose, SWAP
  the OLD owner's object for a mutant, require the answer UNMOVED.** `[M]`
  P4.9b pre-carve it MOVES on every geometry (slab scheme-swap rel **5.000e-02**;
  cyl fp(4,6)/fp(4,8) and sphere closure-swap **4.60e-02 / 5.31e-02 /
  1.196e-01**), so the §6c red-first reading is measured, not argued. ⚠ THREE
  traps, each making it silently green for the wrong reason: (a) mutating ONE of
  the consumed surfaces certifies one route — `[M]` `cell_contribution` alone
  reads `array_equal=True` on all 3 curvilinear rows because `.solve` consumes
  `advance_psi_half` + the minted constants; (b) a DRIVER that re-poses
  internally (`sweep_once` builds the operator at `:814`, i.e. AFTER the swap)
  measures nothing; (c) a surviving cache MASKS the swap, so the gate needs an
  ACTIVATION leg (a freshly posed operator over the mutant hub MUST move).
  → `L64a`
- **⛔ Two of the surfaces a walk "reads off the mesh" can be base
  STATICMETHODS — re-plumbing them is value-inert BY CONSTRUCTION and only a
  structural read-set gate can witness them.** `[M]` `source_emission` /
  `cell_average` live on `DiscretizationSchemeBase`, so `mesh.scheme.X` and
  `op.spatial_closure.X` resolve to the SAME function object. Mode 12 at the
  dispatch. Resolve every mutation through the MRO too — `[M]` 4 of 9 surfaces
  are off the concrete class (2 staticmethods, 1 base function, 1 base
  property), so a concrete-class battery binds **5 of 9** and reports a
  confident partial zero. → `L64d`
- **⛔ A design memo's HAZARD PROSE is a claim — run it.** `[M]` the ruled
  "silent, plausible-wrong k" of a wrong-family closure is LOUD on every
  geometry that matters (sphere `TypeError` naming the requirement; cylinder
  `IndexError`; slab bit-identically inert at `0.0000e+00`, being the default).
  The no-guard RULING survives; its justification sentence must not reach the
  ctor docstring, or it reads as licence to add the forbidden guard. A
  characterization test freezing such a ruling asserts CONSTRUCTIBILITY only —
  one positive leg, no negative — and says why in its docstring (`vv` #11).
  → `L64g`

## 2. Harness discipline — the instrument lies before the code does

**RULE: an all-blind mutation verdict is a broken harness until a positive
control proves otherwise** — `vv` anti-#17 carries the rule and both ORPHEUS
instances (a privately-loaded test-module copy; a summary parser defeated by ANSI
codes). → `L34d`, `L35l`

- **⛔⛔ NEVER quote a NEXUS-derived per-node test count as coverage — in THIS
  codebase it measures the RESOLVER, not the suite.** `[M]` static `callers`
  misses **217 of 229 (94.8 %)** of the tests that actually execute
  `OperatorSum`; the only `OperatorSum` node with any static in-edge is the
  CLASS (39 constructor calls) while `.apply`/`.inverse`/`.apply_transpose`/
  `.assemble` have **ZERO** and all nine members fire at runtime. Tree-wide,
  `[M]` **21.3 % of `calls` edges point into `unresolved`** and the top named
  targets are `op.apply` 265 / `A.apply` 114 / `L.apply` 100. ⭐ The perverse
  part: this is nexus #16's dispatch gap, and its severity is a CONSEQUENCE of
  Cardinal Rule 2 — `coding-elegance` Pattern 1 spells every operation as a
  dunder on a domain type, which is exactly what the resolver cannot follow, so
  **the better the architecture gets the blinder the call graph becomes.** The
  number I nearly shipped: "`[M]` 80.3 % of production nodes have ZERO tests",
  queue topped by `SNSolver.solve_fixed_source` (in-calls **0**, while the free
  function `solve_sn_fixed_source` has **233**). Repair = the runtime overlay
  (`[M]` 163 → 404 nodes reached, `DiscreteMeasure.quotient` 0 → 9 test-callers);
  the only tell that caught it was IMPLAUSIBILITY, not the instrument. → `L55a`,
  `L55b`
- **⛔ A two-sided JOIN needs a denominator assertion on BOTH sides — the
  unasserted side is the one that fails.** `[M]` my pytest↔graph join reported
  a clean, confident `JOIN RATE = 0.0%` because `node_attrs.value` is
  **JSON-encoded** (`file_path` arrives WITH quote characters, so `relpath`
  yields garbage). I had asserted the collected count and not the graph count.
  Caught only because I knew a node existed there. Corrected: 100 % both ways.
  → `L55e`
- **⛔ A wrapping CENSUS plugin must rebind EVERY module that holds the symbol,
  not just the defining one — and its `(*args, **kwargs)` wrapper can BREAK a
  committed gate.** `[M]` #340 N4.7: wrapping `orpheus.diffusion.solver
  .solve_diffusion_1d` reported **`0 entry calls`** on a suite with 34,
  because the tests import from the PACKAGE re-export; the fix is to rebind
  every `sys.modules` entry whose attribute `is` the original (6 bindings),
  asserted at configure time. Separately the wrapper defeats
  `inspect.signature`, so any reference helper built from it (here
  `reachable_knobs`) reds for INSTRUMENT reasons — read that file's colour
  only from a plugin-free run. Validate the decoder by cross-checking the
  census's distinct node-ids against the independent `-W error` red count
  (`[M]` 24 == 24). → `L46e`
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
- **⛔ A numeric table in MY OWN plan is a `[M]` claim — the obvious continuation
  of an integer sequence is NOT a measurement.** I shipped orbit counts
  `{14:6, 16:7, 18:8}`, extrapolated from a committed `1,1,2,3,4,5`; `[M]` the
  truth is `{14:7, 16:8, 18:10}` — the count is `p₃(N/2−1)`, and the sequence
  looks linear for exactly as long as the committed part shows it. Compute the
  extension in the same probe that produced the committed rows, or mark the row
  a placeholder. → `L42d`
- **A pre-existing red must be CHARACTERISED, not counted — otherwise it masks
  the change's own reds.** `[M]` the one red in the very file #337 re-baselines
  is the **GL** case (1 ULP, 8/60 elements) while both **LS** cases pass —
  i.e. it sits on the one rule the change does not touch. That detail yields an
  instruction (`-k "LS4 or LS6"` when attributing; do NOT absorb it into the
  re-baseline, which would hide an unrelated regression inside a legitimate
  one); "1 pre-existing red" yields nothing. → `L42g`
- **⛔ Two overlapping grep predicates in a brief are NOT two work items — compute
  the UNION and print both set differences.** `[M]` CS3's "7 files pinning the
  flux+flux TypeError" + "~16 affine raise-sites" measured as 16 sites in **10**
  files and **8** files by the flux+flux predicate; union **12**, minus 1 false
  positive ⟹ **11** — 4 raise-site-only, 2 flux+flux-only, so either list alone
  misses 2–4 files. And triage a concept grep by MEANING first: `affine` names
  three unrelated things here (the flux torsor; the DD recurrence's affine-in-
  (b, ψ₀) structure; the affine BOUNDARY law `affine-bc-form`, which alone has
  ~18 `:eq:` citers and must NOT be touched). → `L58e`
- **A test-design dispatch's prose grep is a free audit of the PREVIOUS
  change's retirement pass — grep the CONCEPT, then sort hits by TENSE and read
  the JUSTIFYING sentence, not just the number in it.** `[M]` grepping
  `doe=3`/"degree 3" for #337 found 10 stale test comments AND two production
  docstrings already present-tense-false from #327 two days earlier — plus one
  hit that is a load-bearing ARGUMENT ("doe=3 EXACTLY integrates the degree-2
  Y₁·Y₁ moment"), which survives but must be re-stated, never find-replaced.
  → `L42f`
- **MEASURE the reachable subset before rationing the battery — a brief's cost
  estimate names a DIRECTORY, not the tests the carve can reach.** A "≈5.5 min"
  budget measured **9.40 s** for the subset that mattered (35×), turning a
  16-mutation battery from "ration it" into ≈6 min. An over-stated cost silently
  shrinks the battery, which is the same loss as a blind gate and harder to see.
  Also re-measure the pre-declared RED: if the brief names one known failure,
  confirm it is the ONLY one, so every later red is attributable with no triage.
  ⭐ Second measurement, same shape, 2026-08-20: `[M]` `tests/numerics` whole =
  **329.66 s** vs **2.72 s** for the four files a carve could reach (122×), and
  `tests/sn` whole ≈ **80 min** (extrapolated from ~6 % in ~5 min) ⟹ it belongs to
  the pre-merge ≥90-min gate, never to a per-arm battery. The scope that survives
  is `[M]` **1258 passed / 25.30 s**, so a 34-arm battery is ~15 min instead of
  days. **Let the EXCLUDED numbers justify themselves in the plan** — an excluded
  directory with no cost beside it reads as an oversight. → `L36e`, `L60h`
- **⛔⛔ When a behavioural anomaly CONTRADICTS source you have already read,
  re-read it with `inspect.getsource` before believing the anomaly — another
  agent's deliberate mutation is indistinguishable from a production bug.**
  `[M]` I saw `if history.converged:`, then measured a `converged=True` solve
  emitting a warning, and was one step from reporting "the landed commit has a
  live guard bug". The live source read `if history.fully_converged:  # M8
  PROBE — REVERT ME` — the coordinator running MY battery's pre-measurement in
  the shared tree. Reporting someone's mutation as a defect is worse than a
  missed finding. Two habits: **bracket every measurement with the mutated
  thing's state** (`echo GUARD BEFORE …; pytest …; echo GUARD AFTER …` — one
  line, converts an ambiguous number into an attributable one), and treat the
  collision as a GIFT — re-running bracketed handed me M8 for free (`[M]` `130
  passed` flipped vs `130 passed` un-flipped ⟹ **0 reds**, the whole
  justification for the scope gate, measured instead of argued). → `L44h`
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

- **⛔⛔ Read a NULL arm as "my mutation was insufficient" before "the gate is
  blind" — a Pattern-2 TWIN predicate means the survivor guards the gate.**
  `[M]` `SNMesh.reflective_axis_pairs` and `loss_kernel_gauge._reflective_axes`
  are line-for-line the same `len(faces)==2 and all(faces)` test (count vs
  axes); widening ONE → **0 of 25 red**, widening BOTH → **exactly 1**. Grep for
  a second implementation of whatever you just broke. ⭐ Companion: when the
  claim is a CHARACTER identity, compute the stabiliser EMPIRICALLY — the
  canonicality leg `ψ_exact ⊥_G ker A` survives every metric constant across an
  orbit's cells (`[M]` `×2`, `×(1+½ sign Ω·n)`, the partner face's sign, one
  face `×3`: all `≈1.5e-15`; only a random per-DOF metric reds, `1.9e-02`), so
  it is a gate on PARITY, not on metric values — say so instead of claiming a
  sensitivity it lacks. ⭐ And size the POSITIVE CONTROL's regime: a diamond
  weight `w < ½` AMPLIFIES the face mode (`−(1−w)/w = −1.22`) and the arm ran
  25 min for 18 of 25 rows; `w = 0.55` DAMPS it and the same arm is **46 s /
  13 reds**. → `L49c`, `L49d`, `L49e`
- **Run the teeth harness over your OWN new module before delivering it.** It
  flagged a gate I had just written — an activation guard that was a theorem
  about parity, true for every input, surviving authoring and a green run.
  → `L34d`
- **⛔ Run PYRIGHT over your own new TEST module too — it catches the elegance
  defect, not just the type.** `[M]` my first knob-sweep parametrized over a
  string tag + `if entry == …` chain + a `**kw` splat from an untyped dict:
  **24 errors**, all the splat (pyright guessed the wrong parameter). The
  `# type: ignore` reflex would have hidden a real anti-#4 stringly-typed
  dispatch; the principled table-of-typed-callables took 24 → 0 AND let the
  reference read the signature of the very callable the row invoked. Also
  measure the COMMITTED file: mine showed 10 errors, 9 cleared for free by the
  carve, the 10th a **mis-placed `# type: ignore`** (suppression on the call
  line, error reported on the argument line). → `L44k`
- **⛔ Build a source-mutant by TRANSFORMING `inspect.getsource`, never by
  hand-copying the function** — a hand copy is a twin path that drifts from
  production, and a `str.replace` whose target is ABSENT can `raise`, which
  makes the instrument assert its own installation instead of printing a
  banner nobody reads. Smoke-test each mutant's OUTPUT before the battery:
  `[M]` all six printed exactly their intended defect (the partial-carve
  mutant reproduced the welded lie verbatim), so every later red was
  attributable with zero triage. → `L44i`
- **⛔ Attribute an out-of-scope red by AUDITING THE DIFF for arithmetic, not
  by re-running before/after** (which this tree forbids — uncommitted state,
  no `git checkout`). `git diff -U0 orpheus/ | <strip comments> | grep -E
  '=|def |raise |return '` gives the complete added/removed CODE line list; if
  none touches a flux/matrix/residual the change is provably numerics-neutral
  and the red is pre-existing. `[M]` 3 golden-sha reds in a wider run,
  attributed in one command to a Signature-10 stale snapshot on the
  `level_symmetric` family #337 re-seeded — and NOT re-baselined from inside
  an unrelated campaign, which would hide a real regression inside a
  legitimate one. → `L44l`
- **⛔ IMPORT-CHECK the campaign's OWN mutation harness before planning a battery
  on it — a campaign that retires symbols breaks its instruments by MODULE-SCOPE
  BINDING, silently, because nobody runs them between carves.** `[M]` a tracked
  `scratch/mutate_*.py` could not import at HEAD: line 24 bound a private
  predicate the campaign had retired 1 day earlier. Two of its five mutations
  were also stale into anti-#18 (they constructed a keyword a sibling gate
  asserts cannot exist, and declared a pre-carve signature ⟹ their reds are
  `TypeError` crashes, not property reds). One `hasattr` is the whole check;
  repairing the harness is part of the commit, since no negative verdict is
  trustworthy until its control passes. → `L44d`
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
- **⛔ "Prove it does not allocate" cannot be gated by asking a densifier to
  `MemoryError` — `[M]` a 550 GB `np.multiply.outer` is OOM-KILLED (exit 137),
  which fails the RUN, not the test.** Size the gate for SEPARATION instead:
  `[M]` `(2000,)⊗(2000,)` gives dense `32 000 000 B` vs per-axis `32 000 B`
  (**1000×**) in **4 ms**, asserted on reachable `ndarray.nbytes`. Add an EXACT
  structural leg (the dense slot is `None`; no reachable array has
  `size == prod(shape)`) and a BEHAVIOURAL leg — "never densify" implemented by
  DROPPING the metric passes the first two. Rejected `tracemalloc`: the NumPy
  allocator domain makes it version-fragile where `nbytes` is exact and free.
  → `L59c`

- **⭐⭐ Adjudicate a proposed CONSTRUCTION GUARD by INSTALLING it as a plugin and
  counting reds — a per-INSTANCE census (vv#29), never a static site count.**
  `[M]` CS4a, sub-scope baseline **845 passed / 23.89 s**: the charter's
  `space.shape[0] != kernel.ng ⟹ raise` gives **182 failed / 68 errors — 250 of
  845 rows destroyed (29.6 %)**, i.e. unrunnable, not weak; the axis-keyed
  alternative gives **845 passed** with census
  `{checked: 192, skipped_axesless: 578, skipped_nospace: 252, raised: 0}` —
  live on **192 of 1022 constructions (18.8 %)**, inert on **81.2 %**, and it
  **raised 0 times**, so it has no witness anywhere in the suite. ⭐ The
  site-level fraction the design records carried (4 of 13 bindings) understates
  the inertness **4×**: a site census counts call LINES, a running suite counts
  what those lines EXECUTE — put the instance-tier number in the guard's
  docstring. Each plugin asserts its own installation (`RuntimeError` unless it
  binds 4 of 4 classes) and prints a `sessionfinish` census so the DECODER is
  visible. → `L61b`

- **⛔⛔ "Measured-cheap; time it at execution" is an UNMEASURED cost claim —
  price it at PLAN time, and gate it as a COUNT.** `[M]` P4.9b's ruled
  operator-held table: `StreamingOperator` is built **6** times per slab
  eigenvalue solve (**10** sphere), independent of `nx`/`ng`/inner solver, while
  the table is built **exactly 1** time today ⟹ a per-operator memo costs
  `6 × 8.78 ms` on a `284.8 ms` solve = **16.8 %** (24.65 % at 8). Pin the
  builder's call count per solve with the ruled number in the message; a wall
  clock is a flaky proxy for the same question (L24/L25). → `L64b`
- **⭐ Compare an M1 superset PER ARM, never as a union — and read the
  never-red row as an arm-composition gap before calling it blind.** `[M]`
  P4.9b's frozen corpus (27 tests): `m1_scheme` **20** reds, `m1_closure`
  **16**, overlap 10 / 10 scheme-only / 6 closure-only, union 26. The two arms
  PARTITION the corpus by geometry, so a union comparison would hide a
  scheme-side regression behind closure-side reds. The single survivor was 2-D
  wavefront, whose surface (`cell_kernel_batch`) the arm omitted. → `L64e`

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
- **⛔ A witness a brief NAMES for a NEW predicate must be checked at the TIER the
  predicate takes — a cell-level demonstration and a field-level predicate are
  different objects, and the gap is invisible because both are "the negative-flux
  witness".** `[M]` CS3: the cited `TestPositivityFailure` asserts on
  `strat.update(...).outgoing_spatial_flux`, a bare ndarray from ONE cell visit —
  nothing an element predicate on a `Field` can consume. Three probes found a
  strictly better pair through the PUBLIC entry: a converged
  `solve_sn_fixed_source` with `min ψ = −6.40e-01` (2 of 8 entries negative) and
  its benign sibling `min ψ = +2.18e-01`, ONE parameter (`nx`) apart. → `L58f`
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
- **⛔⛔ The OBVIOUS keystone (a derived table) can be IDENTICAL before and after
  the change at a subset of parameter rows — build BOTH configurations and
  tabulate the derived quantity PER ROW before nominating it.** The #337 seed
  change's headline is "the achieved degree rises", yet `[M]` at S12/S16/S18 the
  old and new families reach the SAME degree (11/11, 15/15, 17/17), so the
  degree table is blind to the whole change at 3 of 8 orders — and those are
  exactly the orders where the value gate's tolerance is loosest. A closeout
  saying "the table gates this at all eight orders" is false in a way that
  survives review because the table *looks* complete (`vv` anti-#20). → `L42c`
- **⛔ A brief's REFUSED-family enumeration is a per-family CLAIM to measure.**
  `[M]` "full NODE_ALIGNED product (edge-node + degenerate)" is
  `on_edge_node=True, degenerate=False` — ONE fact. Transcribed, that pins a
  FALSE reason on every product row (L31's blanket-reason trap). Measuring the
  fact-pair per family is what makes "assert WHICH fact fired" writable AND
  makes the two per-conjunct mutations (predicate reads only one conjunct) the
  pair that proves both load-bearing — `[M]` exactly half the refused families
  red under each. → `L43h`
- **⛔ An ISOTROPIC-ansatz row is a provable NON-CATCHER for anything in the ξ /
  azimuthal channel — parametrizing over both builders is necessary, not
  sufficient.** `[M]` the isotropic cylindrical MMS ladder is identical to 3
  s.f. across `folded_product(4,8)/(4,16)/(4,32)`, `product(4,8)` and the full
  staggered rule: the fold is INVISIBLE to it (its source has no `ξ` term at
  all). Declare per ROW which mutation each catches (`[M]` here: iso catches the
  `η→ξ` swap at rel 1.46; only aniso catches `ξ²→ξ` at rel 0.36). → `L43i`
- **⛔ A convergence+value pair can be blind in a band where the error DECREASES.**
  `[M]` scaling the declared `q` by `(1+ε)`: at `ε ≤ 3e-4` the perturbation
  partially CANCELS the `O(h²)` truncation, so `L2(80)` drops to `0.6–0.8×` the
  honest value and BOTH `orders > 1.9` AND any `rtol` value row stay green; the
  order gate first reds at `ε = 5e-4`. So "value + rate" is not a floor — name
  the band and put the real keystone at the tier the defect lives on. → `L40d`
- **⛔⛔ Compute the SYMMETRY GROUP OF THE WHOLE GATE SET, not of one gate — a
  range + an involution identity can be jointly blind to the one flip the seam
  is exposed to.** `[M]` SN's τ was gated by exactly three properties, and all
  three are invariant under `τ→1−τ`: membership `[0,1]` (symmetric about ½), the
  fold box `[¼,¾]` (symmetric), and the reversal identity `τ_m+τ_{M−1−m}=1`
  (`(1−τ_m)+(1−τ_{M−1−m})=1`). `τ→1−τ` IS the march-orientation flip (measure
  the barycentric coordinate from the downstream edge — a one-token index
  drift), and it reddened **0 of those 4 rows**, 6 of 298 tree-wide. Design-time
  and free: intersect the gate SET's stabiliser with the threat model before
  mutating. ⭐ The catcher is a SIGNED law, bit-exact with an exact equality
  case: `(τ_m−½)·μ_m ≥ 0`, `min` exactly `0.0`/`−0.0` at odd N/M (a node at
  μ=0) so `>= 0.0` needs NO tolerance, plus an ACTIVATION leg (`max > 0`, else
  `τ≡½` passes vacuously). ⚠ Never spell it `np.sign(τ−½) == np.sign(cot ω)`:
  `np.cos(np.pi/2) = 6.12e-17 > 0` while `τ−½` is `0.0` there ⟹ FALSE at odd M,
  TRUE at even — a parity artefact reading as a real disagreement. → `L47a`
- **⛔ A parametrize ARGUMENT LIST runs at COLLECTION, so a raising mutation
  reports `FAILED=0` — the flattering direction.** `[M]` building
  `(label, cosines, taus)` in the `parametrize` list called production at import;
  6 of 13 mutations (INCLUDING the positive control) then died as
  `Interrupted: 1 error during collection`, `rc=2`, 2 s, read off the summary as
  "0 caught". Two fixes, both cheap: parametrize by a LABEL and build in the
  BODY; and put `--continue-on-collection-errors` in the battery while counting
  `^ERROR` separately from `^FAILED`. Same family as `L41e`. → `L47d`
- **⛔ A palindromic-index mutation on a symmetric rule is a PROVABLE
  non-catcher — report it, do not gate it.** `[M]` reversing the sphere's
  cumulative-weight order (`w[n] → w[N−1−n]`) is bit-identical (GL weights are
  palindromic): **0 of 298** rows. Sibling of `L43e`. → `L47e`
- **⛔⛔ ASK WHAT FIELD MAKES THE SUT'S OWN RESIDUAL ZERO — a fixture in the
  SUT's kernel cannot rank it, however rich it looks.** `[M]` the shipped
  curvilinear aniso MMS is `A(r) + B(r)η` — affine in the radial cosine at
  every `r` — and the M-M closure is EXACT on `span{1,μ}` BY DEFINITION of τ
  (`4.4e-16` cyl, `8.9e-16` sph, every order; `0` exactly for isotropic). So
  the flagship angular fixture has **zero closure residual** for the very
  scheme it is used to grade. Design-time, one line of algebra, no run.
  Same check kills #319's diffusion-limit instrument for ANGULAR claims (the
  diffusion limit's angular content IS `span{1,μ}`) while leaving it sound
  for SPATIAL ones. → `L48a`
- **⛔ The strengthening axis for an angular-CLOSURE fixture is PARITY, not
  frequency — adding harmonics can LOWER the resolution.** `[M]` reversal
  resolution `1.00×` (m=1 only, BLIND) → `1.33×` (m=1+2) → `1.05×`
  (m=1+2+3, diluted back); `m=3`-alone is blind again. One EVEN harmonic,
  then stop (amplitude saturates at `c₂≈2`). L40b in a new form: the
  τ-independent floor grows faster than the signal. → `L48d`
- **⛔⛔ ONE campaign can carry TWO independent parity rules — name which one a
  fixture is chosen for.** #344's known rule is about the SOLVE (an even first
  axis + a symmetric source leaves the kernel unexcited). The second is about
  the FUNCTIONAL and bites on a pure kernel mode: a face-SUMMED mirror-ODD
  quantity runs the transverse cells against the checkerboard `(−1)^Σi`, which
  cancels when that count is EVEN. `[M]` max face tangential current on a
  unit-norm mode: `(2,2)` **2.5650e-15** (INERT), `(2,3)` `1.9134e-01`, `(3,2)`
  `3.2691e-01`, `(3,3)` `1.1985e-01`, `(3,4)` `2.4252e-01`. The obvious
  "smallest box" is the one that kills the witness. → `L49b`
- **⛔ Re-measure a brief's "the fixture stopped discriminating" number with
  the OTHER axis refined.** `[M]` at `nx=80` the 12-scheme spread is `2.7×`
  and at `nx=320` it is `9.2×` (`n_φ=32`); `1.8×` vs `18×` at `n_φ=64` —
  half the reported blindness was a spatial floor the source plan had itself
  measured two sections earlier. The ladder is also a COST win: the right
  functional was clean at `nx=80`, 8× cheaper. → `L48e`
- **⭐ When a production branch discriminates on a FIELD, CENSUS that field across
  the fixture corpus before writing the gate — a corpus uniform in it makes one arm
  witness-less, and `vv` #17's granularity trap then fires at the FIXTURE tier
  (the arm's mutation has an EMPTY red set and no obvious diagnosis).** `[M]`
  2026-08-20: `get_mixture(region, ng).eg is None` for **all 12** shipped
  `{A,B,C,D}×{1g,2g,4g}` pairs, and for `get("homo_2eg_n2n")` too — so the CS1
  `EnergyAxis.from_grid` arm had no witness anywhere, while `synthetic(ng)` had
  twelve. The witness existed one line away and only a grep for the FIELD (not the
  arm) found it: `tests/homogeneous/test_homogeneous.py:415-417` builds the repo's
  ONLY `eg`-bearing homogeneous mixture, via `dataclasses.replace(base, eg=…)`.
  → `L59b`
- **⛔ Widen the witness census from the CORPUS tier to the MEMBER tier: an
  invariant arm can be unspellable on one member of a union.** `L59b` says a corpus
  uniform in the discriminating field leaves an arm witness-less; the sharper form is
  that a *member* can make it unspellable by construction. `[M]` CS1.5: the
  eg-coherence invariant ("assigned materials must agree on their energy grid") can
  NEVER fire on the infinite-medium member, which has exactly one assigned material —
  its only witness is a structured medium with ≥2 regions, i.e. the arm that ships
  production-unreached. A whole-invariant mutation then reddens via a DIFFERENT arm
  and the run reports "gated". ⟹ for every arm, name the member AND the input;
  an arm with no witness on any member is deleted or declared unfalsifiable in its own
  docstring. → `L60e`

- **⭐ Re-run the P4.9a activation question PER PHASE and PER CLAIM — the answer
  can INVERT inside one campaign.** P4.9a found every frozen artifact blind (a
  congruence-class gate) and the reflex is to carry that forward. `[M]` P4.9b is
  the opposite: `StreamingOperator.__init__` fires in **26 of 27** frozen
  artifacts (so the corpus pins step 1 universally and nothing new must be
  BUILT), while the scan-cache re-pose is activated by only **15 of 27**. Two
  denominators inside one phase; and the two halves of the same step have
  DISJOINT activating configs (`[M]` per-cell scheme dispatch: slab **80**,
  curvilinear **0**; closure dispatch: slab **0**, curvilinear
  **3 192–14 496**), so neither geometry family alone is an acceptance set.
  → `L64f`
- **⛔ A "the reads re-plumb" done-when is DESIGNED-RED when the ruling puts a
  third of them out of scope — partition by ATTRIBUTE and ship the partition as
  an executable read-set gate.** `[M]` P4.9b: `spatial_basis_per_axis` (15) +
  `is_multi_moment` (6) are SPACE facts the hub owns by the same ruling that
  motivates the phase, while the per-cell kernel surface is 9 reads. Wrap the
  hub's objects in a recording descriptor after the pose, run one sweep + one
  matvec, assert the recorded attribute set ⊆ a declared allowlist. → `L64c`
- **⛔ A retirement's MEMO/SLOT inventory is a claim to count, and the contract's
  only witness usually dies with the slots.** `[M]` P4.9b: the design named ONE
  mesh-attr memo; the tree carries **three** (`_geom_cache`, `_coll_cache`,
  `_pole_mirror_cache`), and `_coll_cache` is re-stamped by
  `rebind_cross_sections` — the only reason a σ rebind is not stale
  (`_ensure_coll_cache` never validates σ). The rebind contract's sole witness
  asserts on the retiring solver slots, so its re-pose owes a THIRD leg that
  does not exist today: **net-new teeth created by the retirement itself**.
  → `L64h`

## 4. Reference, claim layer, and the proactive refutation

- **⭐ A test's CLAIM KIND is the PROVENANCE of its expected value — THEOREM /
  REFERENCE / RECORD — and it is a different axis from `l0`–`l3` (which grades
  how GOOD the reference is).** THEOREM = entailed by a law holding for every
  admissible input (identity, adjointness, involution, conservation,
  `M − N ≡ A`); red ⟹ the object violates its own definition and every other
  claim on that subject is VOID. REFERENCE = a structurally-independent external
  route (`vv`'s three pillars); red ⟹ implementation disagrees with the math
  *here*. RECORD = whatever the code produced on a chosen day; red ⟹ *something
  changed*, **zero** information about which side is right
  (`numerical-bug-signatures` Sig 10). ⛔ It cannot be DERIVED — `[M]`
  `assert_allclose` appears in **218 files** and is used identically for
  closed-form, MMS and frozen baselines. ⭐ The audit it unlocks: **every RECORD
  subject must also carry a THEOREM or REFERENCE test** (the hunt task #51 ran
  by hand twice). Honest limit, ship it inside the audit's output: it finds
  subjects with NO independent pin, never a BLIND one — that is mutation's job.
  ⛔ `vv` L51's CONSTRAINT/RANKER/DIAGNOSTIC is the WRONG vocabulary here: every
  collected pytest test is a CONSTRAINT by construction, so that partition has
  one non-empty cell. → `L55i`
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
- **⛔⛔ A brief's headline NUMBER carries an unstated REGIME — reproduce it
  before designing to it, and if it only reproduces off the production path,
  say so.** `[M]` "the angular recurrence produces `min ψ̂ ≈ −77`" reproduces at
  **−76.9** only with a RANDOM ψ and a ZERO seed; on the production value path
  (converged flux + the marched ψ½ STATE) the same fixture gives **+0.1337 /
  +0.1286 / +0.1287** at n_φ = 6/8/16 — strictly positive, within 12 % of
  `min ψ`. The sign is a property of the SEED's consistency, not of the scheme
  (the recurrence's fixed point at flat ψ is `ψ̂ = ψ`). ⭐ Then pin the
  MECHANISM, not the observation: `A(M) = max_m ∏_{k≤m}(1−τ_k)/τ_k` is
  solve-free, a pure function of the chart, has a closed-form independent
  reference, and explains BOTH regimes (`[M]` 2.41/2.73/3.36/…/9.44 at
  M = 2…32 — the 9.44 is the figure `vv` #24b quotes, previously uncommitted).
  Companion: `∏_all (1−τ)/τ = 1` exactly, since the reversal identity makes the
  numerators the denominators re-ordered — and that leg is INVARIANT under
  `τ→1−τ`, so it must be declared a non-catcher for the flip. → `L47b`, `L47c`
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
- **⛔⛔ A brief that says a relocated computation uses "plain / flat / simple"
  arithmetic has named TWO conventions — enumerate the candidate spellings and
  MEASURE the spread before writing the pin; the spread IS the pin's
  discriminating power and its tolerance.** `[M]` CS3's "SI computes the
  diagnostic trajectory from flat norms": interior-leaf `space.norm` → flat
  `np.linalg.norm` is **2.29e-16** (≤1 ULP), but the whole-composite
  `_l2_norm(displacement)` — the spelling the loop already has in hand, which
  additionally ravels the boundary block — is **4.71e-3**. The dangerous rival is
  the CONVENIENT one (it deletes a helper). ⟹ `rtol=1e-12`: 4 orders above the
  harmless rival, 9 below the dangerous one, both measured. ⚠ And a "harmless"
  rival can be phase-scoped: the ULP agreement holds only because
  `inner_product_weights is None` TODAY; under a physical `V_cell·w_n` metric the
  same swap moves ρ by **1.12e-3**, so a LATER phase of the same campaign will
  legitimately red the pin — make that a blocking ruling, not a surprise.
  → `L58a`
- **⭐ A "the fall-through lands on X after we delete Y" claim is a RUNNABLE
  experiment, pre-carve: invoke the base implementation the carve will expose, on
  the inputs the guard must still refuse.** `[M]` `Field.__add__(a, b)` on
  cross-mesh same-shape fluxes already REFUSES (`across distinct SNMesh
  instances`) on all 4 constructible leaves — so the charter's owed negative
  control became a specified test with its `match=` fragment measured. ⚠ Assert
  the DISCRIMINATOR in-test: `a.space == b.space` is **True** (`FunctionSpace
  .__eq__` is `(name, shape)`) while `a.space is b.space` is False, so only the
  MESH arm refuses and a row omitting that precondition degrades silently.
  → `L58d`
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
  ⭐ **Widened, and it recurred within one dispatch:** the ask itself can LAND
  mid-design. `[M]` HEAD moved twice while I wrote (`143e6e2a` → `ce6607f5`) and
  one of three asks — the σ_y-parity gate — shipped, re-posed exactly as I was
  deriving it. Deliverable collapsed from *design a gate* to *audit the landed
  one*, which produced six residual gaps (a better deliverable). So: `git log` +
  `git status` at the START **and** the END, and run an existence check per
  promised DELIVERABLE, not only per named symbol. A useful side-observation:
  the transient state had the module's SECTION-HEADER prose already claiming the
  new behaviour while the function still returned the old one — no gate in the
  tree could see that, which became residual gap "assert the builder ships what
  its docstring says". → `L43a`, `L43j`
- **⛔⛔ "Solve for X on interval I" is a CONSTRUCTION only if the root is unique
  there — run a sign-change scan and report the root COUNT before accepting the
  design.** `[M]` #337's briefed "root-find over μ₁² ∈ (0, 1/3)" has **TWO**
  roots at 4 of 9 orders, so `brentq` **raises before it starts** (same sign at
  both ends) — the sentence is literally unrunnable there. Two roots ⟹ the plan
  owes a SELECTION RULE with its own two-legged gate: (a) the shipped value is
  the selected root, (b) the discarded root is exhibited and measured bad
  (weights −0.6…−7.7). Leg (b) is what makes the rule a reason instead of a
  coincidence — and note the obvious mutation reds by RAISING, so it attributes
  nothing (`L31`/`L25`); only leg (b) is attributable. Also scan the inner
  solve's TOTALITY over the same interval: `[M]` 31 % of the S20 bracket is
  rank-deficient and the residual is DISCONTINUOUS where the row selection
  changes. → `L42a`
- **⭐ When a plan's blocker is "is this small number real or is it
  arithmetic?", answer it with arbitrary precision instead of escalating —
  one probe closed a blocking user ruling.** A frontier decided by
  `+1.75e-4` vs `−2.19e-4` in an ill-conditioned solve looks exactly like a
  float64 artifact; re-doing the WHOLE construction in mpmath at 50/60 dps
  (exact targets, exact Gram–Schmidt, exact LU) agreed to every printed digit.
  Keep the reasoning ✅ ANSWERED, not deleted (`plan-authoring` §3), and pin
  the MARGIN VALUE beside the sign so a conditioning regression reds first.
  → `L42e`
- **⛔⛔ AN ADJUDICATING INSTRUMENT (one that RANKS designs) is a different
  object from a gate and owes four checks, none of which is a mutation:
  BASIS, RANK-CORRELATION, cost-against-alternatives (`vv` anti-#24 a/b/c) —
  and ⭐⭐ the ZERO-SET check: solve `instrument(candidate) = 0` for the
  candidate; if the answer IS the incumbent, it measures
  distance-to-the-incumbent and always confirms what is shipped.** `[M]` #235:
  three live proposals shared one zero set (the diffusion-limit test, the
  shipped MMS ansatz, the η-weighted closure defect) because τ is DEFINED as
  the barycentric coordinate; a fourth (R&L `|τ−½|/w`) has zero set `τ≡½`,
  i.e. it is one CANDIDATE wearing a criterion's clothes. Declare each
  instrument CONSTRAINT / RANKER / DIAGNOSTIC in its own docstring — the
  graveyard died of silent promotion. `[M]` 4 of 6 dead instruments die at
  the `<1 s` solve-free pre-flight. → `L48a`
- **⛔⛔ The GRADED FUNCTIONAL is a design choice with its own stabiliser: an
  INTEGRATED one admits signed cancellation and can rank GARBAGE ABOVE
  PRODUCTION.** `[M]` same solves, `n_φ=64`: scalar-flux `L2` ranks two
  garbage τ permutations `1.6×`/`2.0×` **better** than production and is
  blind to a 2 % jitter (`1.04×`); angular-flux `L2` ranks them `17×`/`8×`
  worse and resolves the jitter at `3.96×` — dynamic range `2.1×` vs `40×`.
  Mechanism, exact: the defect is `∝cos(mω_m)` and `Σ w cos(mω)=0` — **the
  identity that makes the manufactured reference closed-form is the identity
  that annihilates the defect in the graded quantity.** Grade the
  UN-INTEGRATED field when one exists. Every gate in
  `tests/sn/verification/mms/` grades the integrated one. → `L48b`
- **⭐ A CONTINUOUS homotopy beats "rank agreement": require MONOTONICITY in
  `w` along `blend(w)=(1−w)A+wB`** — 5 solves, falsified by one non-monotone
  triple, strictly stronger than a discrete ranking. Stratify the ensemble
  (NEAR = a 2 % jitter · MID = the real rival · FAR = garbage) and require
  the threshold on **NEAR∪MID alone** — a ρ over all three is dominated by
  the garbage split (how `D` scored `+0.75` then `+0.06`). The ensemble MUST
  contain the pair inside the stabiliser you fear. → `L48f`
- **⭐⭐ "This comparison is BELOW MY RESOLUTION, and here is the number" is a
  first-class deliverable, not a failure.** `[M]` best fixture + best
  functional resolves garbage `17–40×`, a 2 % jitter `2–4×`, the
  march-orientation reversal `1.25–1.6×` — and **NOT** the two candidates the
  campaign actually argued about (`1.34/1.10/0.95×`, sign flipping with
  `n_φ`). That makes "decide on constraints + the primary source" the sound
  route rather than a fallback. ⛔ Related: **closure-EXACT is not
  accuracy-optimal** (`[M]` 3.2×/7.8× worse), so every closure-residual
  instrument is a DIAGNOSTIC. → `L48c`, `L48g`
- **⭐ The keystone's ORACLE choice decides whether it catches anything — same
  assertion shape, 8 orders of sensitivity apart.** `γ₋ψ == spec.evaluate(...)`
  is self-consistency: under a magnitude mutation BOTH sides move, green for
  every ε (it still catches delivery-COUNT). `γ₋ψ == the manufactured value
  recomputed from the reference OBJECT` reds at `ε ≈ 3e-12`. For any
  "the answer satisfies the declared condition" gate, ask **which side is the
  thing under test** — if the answer is "both", it is not a gate. → `L40e`

- **⛔ A re-pose can INVERT a migration gate's SENSITIVITY partition — the
  inherited `[M]` characterisation dies by being FIXED, not by being refuted.**
  `[M]` CS1's byte gate was characterised (and the row copied into two other
  plans) as *"BLIND to space weights, LOADED on cell volume; `k_inf` blind to
  both"* — reproduced end-to-end: volumes ×2 ⟹ `flux 397.946→198.973`, rates
  double; space weight ×2 ⟹ **bit-identical**. CS4a re-poses the rate from
  `mesh.volume_measure` to `space.inner_product`; `[M]` the VALUES are 0-ULP
  identical (6 of 6 rows), so the gate stays 8/8 — and the sensitivity **swaps
  sides**. The old anti-claim arm becomes a must-RED arm and a brand-new
  must-stay-GREEN arm appears (the *un-wiring* proof) that could not be stated
  before the change. Run both at BOTH HEADs and put the 2×2 in the docstring.
  → `L61c`
- **⛔ A corpus paragraph can carry an honest `[M]` whose LOAD-BEARING half is
  false, because its experiment varies two things at once — and a `.. implements::`
  directive is part of a re-pose's blast radius.** `[M]` a theory page argued
  *"the quotient point's weight is genuinely consumed, not decorative"* on the
  evidence that the rate functional *"contracts against `mesh.volume_measure`"*,
  with `0.225` vs `0.450` — both numbers reproduce, and the experiment changes
  the carrier volume and the space weight TOGETHER (vv#17's granularity trap, at
  the doc tier), while the separated probe shows the space weight is
  bit-identically inert. ⭐ The carve is what makes the claim TRUE, and
  simultaneously makes the mechanism clause present-tense-false. Same page:
  `.. implements:: <label> :by: <symbol>` onto a symbol the path stops CALLING —
  `dead_references` cannot see it (the symbol still exists; the caller changed),
  so the V&V matrix keeps reporting the equation covered by a transcription that
  no longer runs. → `L61d`

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
- **⛔⛔ For a ROOT-FIND gate the tolerance is `noise / slope`, and the slope can
  collapse 4 orders across ONE parameter family — so a single rtol is a false
  red at one end and a dead gate at the other.** `[M]` #337: the float64 root's
  distance from a 50-digit re-solve runs **1.0 ULP at S4 → 40 653 ULP (8.7e-12)
  at S18** as the residual slope collapses `1.03 → 1.7e-4`. The brief's `1e-7`
  was 4 orders too loose; my own draft's `1e-12` was a FALSE RED at three
  orders. Derive it: `Δx ≈ (evaluation noise)/|f'|` predicted within 5×.
  Tabulate **per row**, ×10, decade-rounded; put the arbitrary-precision value
  in the literal; and STATE what the floor leaves ungated (at S16/S18 a
  sub-1e-10 error is invisible to all three gates). → `L42b`
- **Gating a MORE-accurate implementation against the less-accurate one it
  replaces gates it against the error it exists to remove.** "≤1 ULP vs `np.cos`"
  FAILED at 3.75 ulp — but `np.cos(2πp/q)` is not the true value; against
  100-digit mpmath the new code is 0.57 ulp and the legacy 3.72. State the
  criterion against the arbitrary-precision value, and honestly ("within 0.57 ulp
  everywhere", never "always closer"). → `L34c`
- **⛔ A flat `atol` is wrong in BOTH directions at once — derive it from what
  the quantity DIVIDES BY, and note that two quantities in one seam can need two
  different laws.** `[M]` a τ row carrying `atol=1e-13` *and* a docstring saying
  "a row at N ≥ 64 must widen it" measured **2.247e-13** at N=64 — ~450× too
  loose at N=8 AND a false red at the order it predicted. Derived: sphere τ
  divides an `O(ε)` edge discrepancy by the cell width ⟹ `16·ε/w_min` (`[M]`
  ratio 0.00…3.81 over N=4…128); cylinder τ inherits `cot`'s condition near the
  arc ends ⟹ `40·M·ε` (`[M]` 0.17…7.99 over M=2…64); while the PARTITION the
  same τ reads agrees at a flat **≤1.5 ULP independent of M** (both sides are
  one `cos` of an O(1) angle). ⭐ Same for a negative control's floor: the same
  convention gap SHRINKS like `M⁻²` in edge space and GROWS in τ space, so a
  fixed threshold is a false red in one and a dead gate in the other — derive
  `0.4·sinθ·(1−cos(Δω/2))` for the first, keep a constant for the second.
  → `L47f`, `L47h`
- **⛔ "Bit-identical at the degenerate fixture" is usually 1 ULP — asserting
  `array_equal` on it reds your OWN control.** `[M]` at M=2 the ω-midpoint and
  chord partitions differ by `3.1e-17` (because `np.cos(np.pi/2) = 6.12e-17`,
  not 0). Assert "15 orders below the signal", never "the bits match". ⭐ Such a
  blindness CONTROL still earns teeth: it must red when one of the two
  conventions moves (`[M]` 3 mutations) and stay green when production BECOMES
  the other one. → `L47g`
- **⛔ When a guard compares two independently-accumulated floats, the tolerance is a
  MEASUREMENT over the constructible population PER SUB-FAMILY — never a judgement
  about whether the construction "should" be exact.** `[M]` region interfaces
  (cumulative thickness sums) vs mesh edges: **0 ULP** on slab and on every `uniform`
  mesh, **1 ULP** on CYL/SPH `equal-volume` (`4.441e-16`/`8.882e-16`) — 10 of 4902
  random interfaces, all in one sub-family, because that arm goes through a
  `sqrt`/`cbrt` round-trip the `linspace` arm does not. `==` is a *latent* false red:
  green until someone meshes a curvilinear multi-region geometry by equal volume.
  Ship a derived band (`4 × np.spacing(|x|)`) WITH its discrimination margin (the
  nearest wrong edge is one cell away — 13 orders up), and put the arm that PROVES it
  in the battery: set the guard to `==` and exactly the two curvilinear-equal-volume
  acceptance rows red. → `L60b`

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
- **ADMISSION / REFUSAL carve (a constructor starts refusing a family of inputs).**
  Three archetype-specific traps, all measured on Q5.6's cylindrical-quadrature
  flip. **(a) ⭐⭐ If the input type carries a PROVENANCE field, provenance-keyed
  admission is the cheapest wrong guard — kill it with a two-sided pincer, never
  a single positive.** `Quadrature.quotient()` stamps `folded_by=group`, so
  `if quad.folded_by is None: raise` passes every brief-suggested positive. The
  measured pincer: a HAND-ASSEMBLED rule with `folded_by=None` whose arrays are
  `array_equal` to the factory's MUST construct, and a GENUINE quotient of the
  wrong parent (tag present, `on_edge_node` on every level) MUST refuse. Each row
  asserts its own tag precondition in-test or it silently degrades into its
  positive. **(b) the guard usually makes a production BRANCH unreachable — grep
  the predicate's consumers for `if X: continue`.** `[M]` admitting only
  all-carrying rules made the whole #280 direct-seed fold (its
  `NotImplementedError` included) dead on every constructible mesh, and its
  6-test battery lost its SUBJECT rather than needing a fixture swap. **(c) the
  mesh-tier consequence of admission is the single-source proof** — assert
  `carrying == tuple(range(n_levels))` on the admitted family; it is a theorem
  only if the guard and the consumer read one producer. → `L43b`, `L43g`, `L43d`
- **⛔⛔ "The index addresses the right object" gates: a SYMMETRIC generating rule
  makes half the permutation group invisible.** `[M]` on a σ_y-folded product
  every per-level geometric datum — start cosine, α, ΔA/w, per-level η and w —
  is **bit-palindromic** under `p → n−1−p`, because the GL polar nodes are
  ±symmetric. So the obvious coordination functional (the per-level march-ray
  cosine the seed march divides by) is designed-green for the REVERSAL at every
  mesh and order; only the SIGNED axial cosine and the ordinate index sets
  escape the stabiliser. Ship TWO rows, one per side, and declare the blindness —
  the cosine row still earns its place (it catches every non-reversal
  permutation, `[M]` `s_p` 0.5084 ↔ 0.9404). Enumerate the per-index data the
  consumer reads and test each for permutation-invariance BEFORE nominating a
  functional. → `L43e`
- **RE-SOURCING carve (a consumer's facts move from hand-passed arguments to a
  derived object): ⛔⛔ the blast radius is the READS, not the PARAMETERS, and a
  PARTIAL re-point is WORSE than none.** A brief named three facts moving onto
  `record.first_failure`; `[M]` the consumer body read the old source at **five
  more sites** (the binding criterion — which feeds distance, rate, projection
  AND a branch — plus `min_iterations` and three `n_iterations`). Implemented
  literally, the message welds the inner's knob+budget to the OUTER's criterion
  and projection: every number real, every pairing wrong, and it *looks*
  level-correct. **`awk` the consumer's body for every read of the old source
  before designing gates**, ship the correction as a `line | read | today |
  must become` table at the TOP of the plan, and put BOTH partial mutations in
  the battery (re-point A-only; re-point B-only) — one "did it re-point"
  mutation cannot distinguish them. ⛔ Companion trap: **the existing fixtures
  can make the gate an `X == X` theorem by OBJECT IDENTITY** — on a 1-level
  tree `first_failure is record` returns `True`, so 6 of 7 facts are one
  object's attributes read twice and no input can separate them. The only fact
  a leaf still sees is whatever changed SOURCE KIND (here a caller literal → a
  record field). Two-tier fix: synthetic hand-built nested record (`~0 s`,
  exact geometric trajectories ⟹ analytic `rho`/projection, all facts pairwise
  distinct) as the KEYSTONE, plus ONE cheap end-to-end solve for the thing the
  synthetic cannot see (that production NESTS and stamps the child). → `L44a`,
  `L44b`, `L44c`
- **ONTOLOGY OVERTURN (a role TYPE retires; its operations fall through to the
  base): the carve is byte-identical BY CONSTRUCTION — the arithmetic expressions
  are character-identical before and after — so the whole design problem is
  (a) proving that per CONSUMER rather than arguing it, and (b) the gates the
  RETIRING type's consumers silently lose.** Three moves. **(1)** Enumerate
  consumers and find the existing wall: `[M]` 3 of 7 were already covered by an
  escalatable stored-value gate (`L58c`), and the 2 with NO value gate were the
  accelerated and the ADJOINT paths — the adjoint because it takes a typed
  difference OUTSIDE the iteration loop (`KEigenvalue.measure_stopping_criteria`),
  a site a repair of the loop alone would miss. **(2)** The retiring type's guards
  usually have NO negative test (`[M]` `grep DSACorrection tests/` → 11 hits, ZERO
  `pytest.raises`), so the replacement's teeth are NET-NEW — write them in the
  step that rewrites the guard, never a later one. **(3)** ⛔ A step that MOVES a
  capability must ask what ELSE the moved code did: `[M]` the retiring finder did
  two jobs (walk the composite; decide it carries diagnostics) and the natural
  repair kills the walk — which a scar-tissue gate already exists for. Leaving the
  OLD methods alive is a transient Pattern-2 twin whose committed tests then pin
  the DEAD copy (`[M]` free to avoid: 0 and 1 production call sites). → `L58a`,
  `L58g`
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
- **⭐ PRE-carve dispatch AFTER a grounding pass: re-derive the CONSEQUENCE at every
  enumerated site — do NOT re-verify the enumeration.** A grounding pass is good at
  enumeration and skips per-site consequence, because that needs simulating the carve.
  `[M]` CS1's §P closed 10/10 items and still left four: two "stays on the None path"
  sites that actually gain the new default; two `@verifies` MIRROR tests that
  **cannot** migrate (they build the operand BARE, and no-default-derivation is
  ruled) — one of them the line-for-line mirror of the production line being changed,
  so keeping it pins a RETIRED idiom, the same warrant used to delete the campaign's
  own strict-xfails; a shared xfail-marker constant decorating TWO tests whose reason
  string the carve falsifies; and two refusal messages whose PREDICATE stays right
  while the stated REASON dies. Cheapest form: per site write the one sentence
  *"after the carve this site {keeps/gains/loses} X"* — the sentences that will not
  write are the findings. ⚠ Related trap measured the same day: a composition guard
  that SKIPS a `None` operand lets `M⁻¹(bound) @ F(None)` construct happily with
  `domain = None`, so the breakage surfaces one call later. → `L59f`
- **TYPE-ABSENCE / union restructure (a sentinel becomes a typed union; one member
  stops having attributes):** three rules, all measured on CS1.5.
  **(1) ⛔⛔ `getattr(obj, "attr", default)` swallows the `AttributeError` a refusing
  property raises**, so every duck-typed consumer degrades SILENTLY — `[M]` the
  homogeneous `C` would have gone space-anonymous and un-done a whole landed floor
  with no exception anywhere. Grep `getattr(<receiver>, "<name>"` and `hasattr(`
  and read what each hit DECIDES; the attribute set partitions into *may be absent*
  and *must stay legal*, and the second class is decided by the duck-typed readers,
  not by the concept. ⭐ Corollary for the plan: when a promise LEDGER and the design
  BODY disagree about deleting such an attribute, the ledger is the dangerous one —
  it is the thing someone ticks off. **(2) ⛔ `eq=False` may be FORCED by the field
  types, not chosen** — `[M]` a frozen `eq=True` dataclass holding a `dict` of numpy
  dataclasses has a `__hash__` that ALWAYS raises and an `__eq__` that raises for
  equal-but-distinct values; three lines of toy probe settled a two-document style
  argument, and the test-design consequence is that content identity must be gated on
  what the type MINTS, never on its own equality. **(3) ⭐⭐ Say WHICH measure you are
  mutating** — one field can feed two consumers with opposite visibility: `[M]`
  doubling the SPACE weight leaves `k_inf`/`flux`/`sig_prod`/`sig_abs` bit-identical
  (a rank-1 point axis makes the space measure invisible end-to-end, L59a's dual)
  while doubling the CELL VOLUME halves the flux and doubles both rates; `k_inf` is
  blind to both, being a ratio. One arm per consumer; never credit a k-level row for a
  measure claim. → `L60a`, `L60d`, `L60g`

- **KERNEL/DATUM MINT + CONSTRUCTION BINDING (an operator's data becomes a
  first-class frozen type; the space becomes mandatory):** four rules.
  **(a)** ⛔ The "view over the existing data" it is described as may be a
  **writeable alias with production reach** — `[M]` the carrier's per-material
  accessor returns the CACHE OBJECT (`is` True across calls, for the list and
  its elements), `writeable=True`, and `stack[0][0,0] += 1.0` moves the
  assembled operator by exactly `−1.0`. ⟹ the equivalence gate is **bit-identity
  (`array_equal`), NEVER view-identity (`is`)** — an `is` gate asserts the hazard
  as the contract — plus a **non-aliasing** gate (`is not` + `writeable is
  False` + a carrier-mutation-does-not-propagate leg), whose per-arm proof is
  "copy but leave writeable" reddening the flags leg ALONE. **(b)** ⭐ The honest
  REFERENCE is the ORIGINAL sparse/authored source, not the cache: whichever way
  the design goes (absorb the cache, or delegate to it), a cache-vs-cache gate
  goes tautological (`coding-standards`' single-sourcing clause) while a
  sparse-source gate survives. **(c)** ⭐⭐ "NO apply arm deleted" is a
  **BEHAVIOURAL matrix**, not a grep: `[M]` 4 operators × 3 carriers with every
  cell distinct INCLUDING a `TypeError` cell, because a `singledispatchmethod`
  registry-keyset gate is structurally blind to the operators whose arms are
  `isinstance` chains (`[M]` two of five classes have an EMPTY registry) — and
  those are the ones a "tidy-up" hits first. **(d)** ⭐ A brief's *"flip X too if
  free"* is a claim to MEASURE: `[M]` the annotation-face row does NOT flip when
  only a sibling's space goes mandatory, the ledger moves 16 → 14 (write the
  arithmetic in the commit), and one gate **DIES** rather than flipping — the one
  whose whole subject is the now-unconstructible space-less build. ⛔ Delete it;
  repairing it by adding the new argument turns a real pin into `X == X` under an
  authoritative name. Size the flip by READING the regex hits: `[M]` 10
  "space-less" test constructions were **9**, one being a message string.
  → `L61e`, `L61f`, `L61h`

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
- **⭐ An identity derived from a RENDERING is only as injective as the renderer —
  ndarray `repr` TRUNCATES with `...`, so two distinct long weight vectors give the
  SAME name, and no small-toy gate can see it.** Derive from `.tobytes()` through a
  digest, and state the float caveats (`-0.0` vs `+0.0`: different bytes, equal
  values; `nan`: equal bytes, unequal values). ⟹ **an injectivity gate needs at
  least one pair whose SHAPES are identical**, else `shape` carries the
  discrimination and the NAME is never tested. → `L59e`

## 9. Pointers

- **Characterization vs guarantee:** GUARANTEE tests carry `verifies(...)` and
  assert what IS correct; CHARACTERIZATION tests carry NO `verifies(...)` and
  bound a limitation ONE-SIDED (no upper bound, so a future fix keeps them green).
  To pin a floor a fix claims to remove, measure the floor's SCALING with the
  OTHER axis — `err(S32) < err(S16)/2` is falsifiable where "the floor is gone" is
  not. An out-of-scope defect gets a POSITIVE assert-the-defect gate with a loud
  message, NOT a NON-STRICT or IMPERATIVE xfail (the first flips silently to
  `xpass`, the second cannot flip at all — see §1's deferral-retirement rule).
  A `strict=True` MARKER is the exception and the preferred spelling: its XPASS
  is a FAILURE, so it retires itself. → `L5`, `L16`, `L45`
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

## CS4b additions (2026-08-21) — grouped by the families above

**Gates that cannot red.** ⛔ **A guard hoisted to N call sites needs its witness
table MEASURED, and the arm most likely to be blind is the TRANSPOSE / second-arm
TWIN of a witnessed forward arm.** `[M]` 8 of 22 mesh-identity guards redden
nothing over 3936 rows, and two of them are `apply_transpose` / the second
`.solve` arm whose forward siblings ARE witnessed. A whole-guard mutation
certifies both. → `L62a`
⛔ **A defaulted `getattr` in a guard's condition is a coverage claim with a
hidden expiry date** — the day the attribute retires the branch is unreachable and
NOTHING fails. Grep the retiring name inside `getattr(`/`hasattr(` as the
retirement audit's fourth search. Promoted into `vv-principles` #28. → `L62b`

**Config blindness / claim layer.** ⛔ **"Re-point the space" is not plumbing when
today's space has NO metric** — `[M]` `Field.l2` moves **41 %** (ratio 0.5927),
not ULP, the moment family 18 %. Check `space.inner_product_weights is None`
before believing any re-point is neutral. → `L62c`
⛔ **A phase-ordering hazard written INTO a gate is a claim with a shelf life** —
`test_si_diagnostic_trajectory.py` says "CS2 owns re-deriving these numbers";
CS4b's own mechanism moves it earlier, so the sentence AND its `[M]` blindness row
become false. Re-check WHICH phase owns a predicted red whenever an earlier phase
adopts the mechanism. → `L62c`
⛔ **A space-equality identity doctrine grants permissions nobody enumerated** —
`[M]` besides the argued BC-blindness, the SCALAR family becomes QUADRATURE-blind,
so `φ_S4 + φ_S8` becomes legal. Enumerate the blindness table PER FAMILY, from a
fixture that varies one thing at a time. → `L62f`

**Reference & claim layer.** ⭐⭐ **When one concept has two spellings a constant
apart, do not choose — ask which ARROW each one is.** `[M]` the two isotropic
embeds differ by exactly `Σw`, and the larger one IS the metric ADJOINT of the
retraction (`⟨Rψ,φ⟩=⟨ψ,R†φ⟩` at nulp 1.0) while the smaller is the SECTION
(`R∘E=id`, bit-exact). Naming them apart makes the missing-factor class
(ERR-051) unspellable. → `L62d`

**Harness discipline.** ⭐ **Build the attribution scope FROM the positive
control's red set.** `[M]` the all-arms control's 29 reds defined a 548-row /
3.31 s scope, so a 21-arm loop cost 70 s instead of 9 min. → `L62a`
⛔ **The inherited phase battery is scoped to the PREVIOUS phase's blast radius.**
`[M]` 7 of 29 reds and the gate the carve is predicted to break are outside the
14-path CS4a scope; the amendment costs +175 rows / +2.90 s. Re-derive the scope
from THIS phase's cone. → `L62f`
⚠ **Check the collection-kill hazard, do not assume it** — `[M]` 0 module-scope /
parametrize-list constructions in the scope; I had written the hazard as live and
my own measurement refuted the sentence. → `L62f`

**Snapshots & exactness.** ⭐ `[M]` gate-ready bit-exact laws found by probing
rather than assumed: `R∘E = id`, `E∘R` idempotent, `R.H == Σw·E`,
`HarmonicFrame.analyse(ℓ=0) ≡ integrate_angular` — all `np.array_equal`. Probe
the algebra before choosing a tolerance; four of these needed none. → `L62d`
⛔ **A derived space NAME is a future landmine when axis identity is
per-SUBCLASS** — CS2's axis classes will change every `of_axes` name. Pin axis
CONTENT and relative identity, never the name literal. → `L62f`

**Config blindness.** ⛔⛔ **When a code path is gated by a parameter's CONGRUENCE
CLASS, the frozen corpus samples one class only.** `[M]` the whole P4.9a carve
runs on cylinder `n_phi ≡ 2 (mod 4)` and NOWHERE else (0 `DD.update` calls on
slab/sphere/`n_phi≡0`, 13 760 at `fp(4,6)`) — so `4, 8, 16, 32` reads as a
refinement ladder and is a single residue. **Every** frozen artifact was blind,
including the plan's own named canary. Run a counting spy and confirm the changed
line EXECUTES before crediting any snapshot as a carve's anchor. → `L63a`
⭐ And the tree usually already knows: two authored comments stated the rule
exactly and shipped the activating fixture — authored, and never adopted by the
frozen corpus. → `L63a`

**Reference & claim layer.** ⛔⛔ **A "delete the twin" carve must ask which
ARITHMETIC FORM the destination spells, not just which module owns it.** `[M]` a
third live spelling sat outside the done-when's grep scope, algebraically equal
and **204 ULP** away (bit-equal 59 % on real τ; 100 % only where `τ == 0.5`
bitwise, which holds on just 2 of 4 degenerate ordinates at one config). Routing
through it costs 1–2 ULP in keff and breaks `array_equal` on 3 of 4 configs for
zero gain. → `L63b`
⭐ **A gate's BUILDABILITY is a legitimate design constraint** — the charter's
`is`-identity gate is unwritable unless the branch calls a closure METHOD; say so,
because it converges with the arithmetic argument on one ruling. → `L63b`

**Gates that cannot red.** ⛔ **Shedding a protocol field disarms every guard that
KEYS on it.** `[M]` the #158 curvilinear refusal keys on `angular_upstream is not
None` and its only witness constructs that field directly — both die with the
shed. Grep every removed field as a GUARD PREDICATE (`is None` / `getattr(…,
default)`), not only as a read, and land the re-key + its witness in the SAME
commit. → `L63c`

**Harness discipline.** ⭐ **Run the mutation BEFORE writing rewire
prescriptions** — a per-gate claim-class verdict guessed from reading is wrong
both ways. `[M]` 25 of 33 reds were a `TestResidual` family the brief never
listed (the real cross-helper gate), while a LISTED gate red **0** because it
feeds the helper's output to the SUT (`vv` #22 shared input). → `L63d`
⭐ **A monkeypatch-only battery is crash-safe BY CONSTRUCTION** — nothing on disk
to restore, strictly stronger than copy-aside + `diff -q`. → `L63h`
⚠ **zsh does not word-split unquoted `$VAR`**: `pytest $SLICE` collected 0 tests
and read as a clean run. Use a driver script with `"$@"`. → `L63h`

**Snapshots & exactness.** ⭐ **The cheapest anchor is a PARAMETRIZE ROW on the
harness that already has the right regime.** `[M]` +1 row (`fp(4,6)`) on the
existing het-σ_t/random-source snapshot gate = **2.1 ms** and 32 M-M calls. Land
it BEFORE the carve (else the snapshot inherits the new code), and ADD a member —
never retune a shared literal, which silently re-baselines the sibling row. →
`L63f`
⛔ **Retiring a helper breaks DECLARED provenance edges — compute the set
DIFFERENCE.** `[M]` one of four `implements` edges was unique to the dying symbol,
and one equation's ONLY three claiming tests were the twin's catchers. → `L63e`
⭐ For a "hand it the constant" move the realistic defect is the CLEANER algebraic
spelling, `[M]` 1–2 ULP — gate with `array_equal`; any tolerance ≥ 1e-15 is a
non-catcher. → `L63g`

## CS5 additions (2026-08-29) — grouped by the families above

**Gates that cannot red.** ⛔⛔ **When the step-1 code is ALREADY in the tree, run
the battery against the EXISTING suite FIRST** — the question is not "will my new
gates redden" but "what does the tree already catch". `[M]` 3 arms over a
184-test anchor set: **0 genuine catchers**, and the ONE red was a cross-process
determinism gate whose subprocess leg is unmutated (a universal digest tripwire,
not a provenance catcher). Banking it would have inflated coverage by one and
hidden that the true count is zero. → `L65a`
⭐⭐ **When a re-point makes `new_path is old_path`, EVERY value gate is a
tautology** — green before, after, and under a PARTIAL re-point. The keystone must
be a ROUTE gate: install a **decoy that is invisible to identity and visible only
through the accessor** (`[M]` weight-preserving, `np.roll`ed nodes: space name and
`==` IDENTICAL, `mu_x` DIFFERENT), then require the answer to MOVE, with an
anti-dud control and a leg asserting the decoy did NOT move the space name.
⚠ `[M]` `-nodes[::-1]` on Gauss–Legendre is the IDENTITY — print the decoy's
discriminating array before trusting the gate. ⛔ **SHARPENED at the
P4-remainder: the rolled-node recipe is REFUSED by production's own α-dome
admission guard on every curvilinear chart — a decoy must clear the guards of
the arm the gate lives on, not merely discriminate. `nodes × 0.9` is the one
that clears both tiers.** → `L65c`, `L66b`
⭐ **Look for the §6c witness that SHIPS before manufacturing one.** `[M]` three
production sites stay generator-less at landing (homogeneous counting point, every
`EnergyAxis`, the MODAL moment axis), so the refusal guards a live state. Found by
reading the field off real objects, not off the design — and record the witness's
SHELF LIFE (one retires when the MODAL half lands). → `L65f`

**Config blindness.** ⛔⛔ **A mint that consumes a FLAT collection is rank-1 by
construction, and d=1 hides it.** `[M]` `DiscreteMeasure.axis` → `(n_points,)`;
against the shipped spatial axis, d=1 is IDENTICAL and d=2 is
`spatial(12,)#3712…` vs `spatial(3, 4)#1dcb…` — a moved space NAME, i.e. moved
space identity for every operator keyed on it. The congruence is **RANK**, and
every axis-suite fixture is 1-D. ⟹ ask what the CONSUMER's rank is before writing
the gate. → `L65d`
⭐ **A branch added to dodge that (`isinstance(mesh, Mesh1D)`) creates the
P4.9a blindness** — the new path runs on one carrier kind only. Parametrize the
gate over the **BRANCH**, and make the other arm assert the OPPOSITE claim
(`generator is None`, digest unmoved). Two arms, two claims; one-arm gates let the
sibling drift. → `L65d`

**Reference & claim layer.** ⭐⭐ **Before excluding a field from an identity key on
DOCTRINAL grounds, check whether the exclusion is also MANDATORY** — that is the
stronger, more durable gate. `[M]` including a `Quadrature`/`DiscreteMeasure` makes
`__eq__` **RAISE ValueError** and `hash` **RAISE TypeError** (ndarray `eq=True`
dataclass; `frozen=False, eq=True` ⟹ `__hash__ = None`), not merely disagree. Pin
the REASON with `pytest.raises` legs on the generator TYPES. → `L65b`
⭐ **The intrinsic law of a provenance accessor is the SECTION law**
`a.generator.axis(a.label) == a` — `[M]` holds 4/4 angular rules and is exactly
what fails at d≥2, so it is the rank blocker's standing gate. State it over
*generator-ful* axes only, or the shipped generator-less sites read as violations.
→ `L65e`
⚠ **Name which half of a comparison is real.** `[M]` both the angular and spatial
"mint vs literal" gates read the SAME array object on both sides
(`q.weights is q.measure.weights`; `volume_measure.weights is carrier.volumes`),
so they pin THREADING (label, shape spelling, `kind`, wiring) and never the
values. The honest digest gate rebuilds the **pre-change literal space IN THE
TEST**. `[M]` the chain has no independent literal anchor at any rank — the
non-uniform `V` vector lives in a COMMENT. → `L65g`

**Harness discipline.** ⚠ **After adding a field to a type, grep the tests for
REFLECTION walkers** (`vars(`, `asdict`, `fields(`) — `[M]` exactly one here, and
it happened to drop the new field only because it descends into
`ndarray`/`tuple`/`list`/`Axis` and nothing else. A walker over arbitrary objects
would have swept the generator's arrays into a no-densification count and reddened
for an unrelated reason. → `L65h`
⭐ `[M]` `dataclasses.replace` round-trips a frozen canonicalizing dataclass
**bit-identically** (bytes, read-only flag, `eq`/`hash`, kw-only subclass fields,
idempotent canonicalization) — so a `replace`-based field upgrade is Pattern 4∩2
safe and needs no hand-written constructor. → `L65h`
⚠ **A `@property` mints a FRESH object per access** — `axis.generator is
mesh.volume_measure` is a latent false red. Assert type + content. Same trap
`EnergyAxis`'s docstring already records for `Mixture.energy_grid`. → `L65g`

## P4-remainder additions (2026-08-29) — grouped by the families above

**Gates that cannot red.** ⛔⛔ **A decoy must survive the PRODUCTION ADMISSION
GUARDS of the arm the gate lives on** — one level past CS5's "print the decoy's
array". `[M]` the CS5-prescribed weight-preserving `np.roll` decoy is REFUSED by
`angular_redistribution`'s α-dome guard on every curvilinear chart (the dome
closes iff `Σ w_m µ_m = 0`; a roll breaks the (µ,w) pairing) and refused again
at the M-M mint (τ ∉ [0,1] / ω-arc). The one decoy admissible at BOTH tiers is
`nodes × 0.9`: axis `==` preserved, `|µ|` moved 4/4 GL4 · 8/8 GL8 · 8/12
fp(4,6), mint ADMITTED with τ + µ_x MOVED. ⛔ `nodes[::-1]` and `-nodes` move
**0/N on every GL rule** (`abs_mu` takes the modulus, GL nodes antisymmetric),
and `gauss_legendre(2)` moves **0/2** under the roll — which is exactly what
`make_tiny_spherical_sn_mesh()` uses. → `L66b`
⭐ **The COURIER'S REMOVAL, not a gate, is what makes a partial re-point
unspellable** — delete the twin field and "did every read move?" becomes
structural. Promote the `dataclasses.fields` name-set row (set EQUALITY, so a
re-addition also reds) from nicety to keystone-support. → `L66a`
⛔ **A member whose only read is a COUNT is a structural non-catcher for any
same-N decoy** — `[M]` `IdentityAngularClosure` reads `mu_x.size` and nothing
else, so it mints bit-identical constants under every same-N decoy; its
discriminator is a DIFFERENT-N axis. → `L66b`
⛔ **The silent mode is the WRONG LABEL**: an object-held axis makes
`axis("ordinate")` value-identical and reds nothing until a consumer looks it
up by name. Gate it with `op.axis == hub.space.axis(label)` — `==`, never `is`
(`Quadrature.axis` is a METHOD minting fresh: `[M]` `is` False, `==` True).
→ `L66e`

**Harness discipline.** ⛔⛔ **A class-NAME §6b census misses a SUBCLASS that
inherits `__init__` and a `Base.create(**kwargs)` registry call** — `[M]` the
banked 8 closure-ctor sites are **10**, and #9 is the previous phase's own
KEYSTONE mutant factory (`_MutantMM`), whose omission kills that gate at
COLLECTION rather than weakening it. Recovery: transitive base resolution over
every `ClassDef` + `super().__init__` + registry `create`; and a `create`
call makes the new parameter's NAME an API surface. → `L66c`
⭐ **Re-running a prior phase's superset battery pins a DIFFERENT claim here:**
it is the INVARIANCE witness (the red SET must be EQUAL per arm), because a
same-sized-but-disjoint set means the re-source landed on a different instance
— the one failure mode no value gate can see. Carry the recorded arm gap
forward rather than banking 26/27. → `L66c`, `L66j`

**Config blindness.** ⛔⛔ **A second MINT SITE hides on the branch where the
producer does not exist** — `[M]` the d ≥ 2 Cartesian arm builds the closure's
operands with `reduced is None`, so it must mint its own axis: a Pattern-2 twin
one label typo apart, on the exact branch the inherited battery's corpus cannot
redden. Parametrize the coherence gate over the BRANCH; the arms assert
different things. → `L66d`

**Reference & claim layer.** ⚠ **The obvious oracle migration DEMOTES the
comparisons it touches** — re-pointing a reference helper to the accessor
production just adopted makes a wrong mint move both sides. `[M]` free to avoid
here (`quad` already in scope at all 10 call sites); the migration changes no
value today, so the entire gain is surviving as a wrong-mint catcher. → `L66f`
⭐ **A provenance UNION cannot answer the consumer's reads, so the refusal must
NARROW** — `[M]` `Axis.generator: DiscreteMeasure | Basis | Quadrature | None`
and `DiscreteMeasure` has none of `mu_x`/`level_indices`/`eta`/`mu_z`/`N`. That
makes the G5-style refusal load-bearing (it is the narrowing) and forecloses
the `# type: ignore` reflex; the pyright ratchet gates BOTH directions.
→ `L66g`

**Snapshots & exactness / cost.** ⭐ **Price field-vs-property at PLAN time; the
existing COUNT gate is blind to it.** `[M]` `Quadrature.axis()` mint 9.61 µs ×
**320** `streaming_terms` calls per solve = **2.49 %** of a 0.123 s solve, and a
fresh-minting property makes `op.axis is op.axis` FALSE (latent false red).
`CollisionCache._build_count` counts a different cache. Owed row: the
field-identity leg. → `L66h`

**Sequencing.** ⛔ **A route gate whose SUBJECT is created by the step it gates
cannot be a later step** — scheduling it after leaves an ungated interval and
makes its §6c red-before reading permanently untakeable. Fuse them; the
pre-carve evidence that survives is the SIMULATION (substitute a decoy for the
soon-to-die twin field). → `L66k`
