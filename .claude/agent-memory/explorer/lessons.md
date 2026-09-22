# Explorer — Lessons (the hot digest)

Behavioural corrections only: *what mistake did I make exploring, and what
changed how the NEXT exploration is run?* Each entry is an imperative, its
tell, and a pointer. The war story behind every entry is cold, in
`lessons_archive.md` under the same L-number: open a section only when a
pointer needs checking. The HOW of every Nexus tool is the preloaded skills
(`nexus-exploring`, `nexus-guide`); the standing directives are `AGENT.md`
Operating Principles 1–7; neither is restated here. L-numbers are stable
(topic files cite them); an entry that retired into a directive keeps a
one-line stub. Rules reach this agent only through a brief's "Rules that
apply to you" line (`articulation`, `code-search`, `coding-standards`,
`instrument-doctrine`, `process-discipline`); an entry that restates a clause
of one of those rules which the brief line does NOT carry is marked
`[uplift]` and stays until the rule's brief line carries it (none is marked
today: the 2026-09-21 candidates landed).

## The six meta-lessons

- **M-1 Every inherited datum is a claim with a vintage.** An issue body; a
  brief's timeline, type table, exemplar, count or `Class.attr (file:line)`
  citation; a plan section marked "retained"; a docstring saying "the ONE
  site" or "X handles it"; a stored numeric tag; a test's self-description;
  a `[M]` on a NEGATIVE ("no consumers", "discarded", "cannot express"). Each
  is verified by its cheapest decisive probe before anything is built on it,
  and the strongest-looking ones expire first because nobody re-checks them.
  → now in `AGENT.md` OP5 (2026-09-21). Instances: L-002, L-011, L-016,
  L-017, L-020, L-021, L-022, L-023, L-025, L-027, L-030, L-034b, L-035,
  L-036, L-045.
- **M-2 The tree moves while you audit.** Open with `git status --short`,
  `git diff --stat` and `git log --oneline --since=<the SECTION's vintage> --
  <scope>`; close by re-running verbatim every search whose EMPTINESS is a
  finding and `git ls-files --error-unmatch` on every file you call landed;
  read the NEW module an intervening commit added and ask each negative
  claim against it. A "zero consumers / cannot express / not yet built"
  verdict is the most perishable finding there is. The brief template carries
  the open/close protocol since 2026-09-21. Instances: L-007, L-012,
  L-020, L-027, L-035, L-036.
- **M-3 A census is an instrument: split its populations before counting,
  and control its zero.** The brief line carries X1/X2 (predicate, tree,
  exclusions, `k of N`, a positive control per shape). The splits this
  corpus keeps teaching: executable call vs prose citation vs duck-typed
  surrogate; receiver-resolved vs bare token; xref vs promissory claim vs
  convention statement; anchored vs unanchored (the English-word substring,
  the same-package homonym); the AST node-type viewport; two counters at
  different depths so "built but never reached" is a visible row.
  Instances: L-009, L-017, L-027, L-030, L-037, L-038, L-039, L-040, L-043,
  L-045.
- **M-4 The graph is blind at known seams; there grep/AST is PRIMARY
  evidence, not a cross-check.** The skill states the general fact
  (`nexus-exploring` "⛔ The caveat on each answer": properties, dunders,
  callbacks and polymorphic dispatch mint no edge; `callers` empty ≠ dead).
  The seams measured here: dataclass FIELDS (no node at all); METHODS
  (`callers` → `nodes: []` + `unresolved`); `@singledispatchmethod` arms
  behind an `apply = _apply_impl` alias; Protocol-typed receivers; function
  objects captured in a dataclass field or a catalogue dict; `@property`
  bodies and class-body installs; docstring `:func:` roles that mint a
  `references` edge (flattery, the dangerous direction); a module's OWN
  import binding (patch the caller's, not the definer's). Instances: L-001,
  L-009, L-017, L-024, L-038, L-039, L-044.
- **M-5 A behavioural question is answered by a RUN on the discriminating
  input with a control beside it, never by a read.** "Does it break", "what
  is in scope", "is it bit-identical", "is the degraded state deliberate",
  "is it hashable", "does the guard bite": swap it and run; spy the frame
  locals; solve the counterfactual; `hash(a)`, `a == b`; ULP-probe a random
  operand. The control: the free baseline (a random orbit with the same
  declared symmetry), the fixture that BREAKS the property, the `None` arm
  and the same-data rebuild, the production quadrature rather than the slab.
  An all-green run may have measured INERT: find the gate and check the
  path routes through it. → now in `AGENT.md` OP8 (2026-09-21). Instances: L-010, L-013, L-016, L-018, L-021,
  L-024, L-025, L-026, L-028, L-041, L-042, L-044.
- **M-6 A verdict names both arms and the discriminator, and hands the value
  judgment up.** Retire vs keep-as-anchor; presence probe vs width probe;
  fold the algorithm vs fold the state; refuse vs reduce, priced against
  the incumbent; the incumbent's implementation vs the challenger's
  fragments; "conflict to adjudicate" when a plan prescribes what a ratified
  counter-design already replaced. The rejected arm carries its structural
  reason (`process-discipline` "A refuted candidate is first-class output").
  Instances: L-004, L-006, L-015, L-021, L-029, L-034, L-036.

## The lessons

- **L-001** → `AGENT.md` OP4 (the four searches). Retained sharpening: two
  constructs make `callers()` lie systematically — `apply = _apply_impl` with
  `@singledispatchmethod` arms, and Protocol-typed receivers
  (`solver: EigenvalueSolver`). A dispatch arm's liveness is decided by the
  ACTUAL input type at the production call site: read the chain, let grep
  enumerate the `op.apply` sites. → M-4; archive §L-001.
- **L-002** → `AGENT.md` OP5 (verify the issue's premise; stale ⟹
  "CLOSE-VERIFY"). Archive §L-002 holds the worked examples.
- **L-003** → `AGENT.md` OP7 (durable shape leads, `file:line` is
  re-derive-via-Nexus). Home placement: durable shape into `AGENT.md`'s
  durable-shape section; a transient line map into a topic file stamped
  with its HEAD, deletable once the campaign merges. Archive §L-003.
- **L-004** A carve verdict is "RETIRE-eligible BY <discriminator> — AND the
  documented KEEP — AND the call is the user's": same math available via the
  surviving helper ⟹ retire; a named independent consumer need, even a
  future one ⟹ keep-as-anchor is defensible. Never pre-decide a retirement
  that turns on "will a future issue consume this". → M-6.
- **L-005** → `AGENT.md` OP6 (git is authoritative for merge status; also
  `process-discipline` "Trust git for merge status"). Archive §L-005.
- **L-006** "Collapse these N shape probes into one predicate": classify each
  as boolean PRESENCE (swap to the typed predicate) or integer WIDTH (an
  honest count feeding allocation — keep), and check whether the probe site
  even HAS the typed object in scope; report "(B) small plumbing", not "(A)
  clean swap", when it does not. → M-6.
- **L-007** Uncommitted edits in the audited subsystem ⟹ re-run the census as
  the LAST step; diff the uncommitted hunks against the brief's items (one
  may be mid-fix, flipping that deliverable to "confirm the in-flight fix");
  stamp line numbers "at final read; tree moving". → M-2.
- **L-008** zsh: an unquoted word starting with `=` (`echo ===`) triggers
  `=cmd` expansion and aborts the WHOLE compound command, silently losing
  every grep sequenced after it. Quote separators (`printf 'NAME\n'`,
  `echo "---"`). → a `code-search` check since 2026-09-21.
- **L-009** A dataclass-FIELD rename is a grep problem: Nexus mints no field
  node (0 of 75 sites came from the graph). Grep
  `\.<field>\b|<field>=|<field>:` plus `replace(`, `getattr`, `asdict`. Before
  any replace strategy, check whether the old token is a substring of an
  English word (`gains` ⊂ `against`: 679 hits; `loss` ⊂ `lossless`, `role` ⊂
  `payroll`): grep with and without `\b` and hand the delta to the
  implementer as a hazard line. → M-3, M-4.
- **L-010** "X is the complement of Y, hence X == Z" needs the split's
  PREDICATE to be exhaustive, and a strict `<`/`>` pair with an epsilon never
  is (tangential ordinates: 4 of 8 on the production cylinder; rank 18 vs 6
  on one face). Re-run across the production data (`Quadrature.product`,
  Lebedev); the slab is the degenerate case for nearly every SN index
  question. → M-5.
- **L-011** A docstring that DELEGATES ("the sweep handles it via its
  face-pair indexing") is the highest-yield falsity shape: grep the invented
  MECHANISM noun, not the symbol (40 healthy hits vs 2, both the claim
  itself); the sibling implementation that REFUSES the same input, with its
  reason, is a free oracle; a strict `xfail` naming the gap outranks prose;
  re-locate quoted prose before judging it (it may be half-stale). → M-1.
- **L-012** On a "blast radius ahead of a carve" brief: `git status --short`
  + `git diff --stat` FIRST and again at close; `git ls-files
  --error-unmatch` on every file you call landed (three "tracked" modules
  were `??`); if the carve is underway, keep the audit as taken, add a
  reconciliation by RUNTIME PROBE against the final tree (a diff hunk left
  `is_adjointable`'s flip ambiguous), and lead with the still-open items.
  Re-run the emptiness greps verbatim at close: a "zero consumers" verdict
  flipped to "consume the rule that landed mid-dispatch". → M-2.
- **L-013** "What breaks if this numeric primitive changes?": swap it (a
  `pytest_configure` plugin) and run the consuming suites — ~200 grep
  candidates became 2 measured items. Classify by "is the RHS FROZEN"
  (`.npz`, a literal, a hash), never "is the comparison exact": same-process
  route-equivalence is immune. Patch every re-export AND every captured
  function object (`object.__setattr__` on the registry field), and print a
  confirmation line. Run the sibling suites too (the only moving snapshot
  was reached through `_generate_snapshots.CASES`). A guard's FIXTURE
  ENUMERATION is where vacuity hides (the eps-gap gate enumerated GL+Lebedev,
  not `product`); check whether the "new" hazard already exists on the
  already-exact sibling (ties: 18–216 per LS rule today). → M-5.
- **L-014** Adjudicating an algorithm against the literature: read the
  paragraph that PRODUCES the equation (Hébert defines `α ≡ 𝒲·η` two lines
  above the recursion — the definition turns "convention?" into a decidable
  closed-form check) and grep the sidecar for the INDEX-DOMAIN prose
  (`interval`, `normalized to`, `for m = 1 … M`). A cumulative recursion
  telescopes under every permutation, so its closure/sum/step gates are
  ordering-blind (`vv-principles` Mode 12). Find the closed form before
  endorsing an MMS run.
- **L-015** "This DOF is redundant, fold it": enumerate the FUNCTIONALS of the
  state (moments, currents, leakage, inner products, the adjoint metric) and
  classify each integrand's parity under the group — odd-parity ones are OUT
  OF THE SPACE, not inaccurate (a fold that was 5e-16 on every even moment
  turned the odd one into +2.94). Split fold-the-ALGORITHM (lift partners
  back; symmetry by construction; no memory win) from fold-the-STATE before
  scoping. When both candidates satisfy the issue's criterion, find the
  structural predicate the tree already keys on (`0 < tau_raw[0] < 1`). → M-6.
- **L-016** A stored numeric tag (exactness, order, rank) is a claim:
  brute-force sweep it, and FIRST measure what a structurally-trivial object
  with the same declared symmetry gives (a random `O_h` orbit was degree 3 —
  identical to the real rule, so the construction contributed nothing).
  Read the WEIGHTS when a moment claim fails (`n_distinct_weights == 1`).
  `assert x.tag == expected` is not a property test: every such line is
  evidence the property is UNTESTED. A `min()` over incommensurable units may
  hide a correct unstated mapping; when a flag's docstring and its registry
  ENTRIES disagree, measure — the entries were right. → M-1, M-5.
- **L-017** Before counting a retirement, a two-number probe: `grep -c
  '<name>'` vs the anchored form (`[^.]name(`, import lines) — a same-package
  homonym (`Quadrature.gauss_legendre` classmethod vs the module function)
  inflated 570 lines to 2 files; hand the anchored pattern to the
  implementer. A test's docstring is not what it pins: resolve the RHS one
  hop (two spellings converging on one implementation are route-equivalence
  however "load-bearing" the docstring says it is); the real frozen
  baselines are `find tests -name '*.npz' -o -name '*.npy'`. A function
  object captured in a dataclass field is a live consumer with zero graph
  edges. → M-3, M-4.
- **L-018** A hand-written lookup table behind a computed fast path has dead
  rows: per row, "which branch answers this?" (2 of 5 live). A bare
  `return False` fallthrough gives a new tag a wrong-but-silent answer. A tag
  with two dispatch branches is invisible on a fixture both accept: find or
  build the input that BREAKS the property (shipped: `product(4,3)` is closed
  under σ_z and not σ_x). → M-5.
- **L-019** Hunting a hidden transformation: grep the chart-defining
  ASSIGNMENTS (`cos_theta =`, `= arctan2(`, `nodes[:, k]`), reconstruct the
  implied frame matrix, and ask the checker whether it can NAME it
  (`cos_theta = mu_x` is a 120° `O_h` element no tag expresses). COUNT a
  partition's parts before believing its name (`octants` → 26 on
  Lebedev(17)). The tolerance-family census is free (three idle epsilons,
  one justified by a symbol that exists nowhere — L-011's shape in a `#:`
  comment). Evaluate a "for the degenerate case X vanishes" docstring ON the
  degenerate case (slab `m>0` slots ≈ 0.83 at ℓ≥2).
- **L-020** The BRIEF's timeline is a claim: `git log -1 --format="%h %ad"
  --date=iso <hash>` per named commit plus the document's mtime (`stat -f
  "%Sm"` on macOS). A "cannot express X" verdict expired in six hours when a
  sibling campaign landed the substrate (`RigidMotion`): read the new module
  and ask the negative claim against it — the ELEMENT exists, the TAG cannot
  name it, the certifier correctly rejects it. N spellings of one concept
  may differ in DOMAIN/CODOMAIN: count the tiers before endorsing a
  unification; the real twin is the boring shared primitive. → M-1, M-2.
- **L-021** A brief's (or page's) typing table `A --f--> B` is a claim about
  MATERIALISED objects: first ask "is `f` an object?" (one grep of the
  accessor's return type) — two of four rows had no construction site, one
  was unfillable in principle. An all-green "does it break?" run may have
  measured INERT: find the gate (`OperatorSum.__init__`) and check the path
  routes through it; report `inert` vs `verified`. An earlier PHASE of your
  own campaign expires like a sibling's. The cross-face threading you would
  "have to add" is often returned by a helper and discarded at the call site.
  → M-1, M-5, M-6.
- **L-022** "What REMAINS of issue X": grep the campaign's own interim
  vocabulary ("today still", "remaining half", "until then", "not yet", the
  `.. caution::` directives in touched modules) and re-verify each hit — the
  honest mid-flight notes are what nobody returns to (four falsified in one
  pass). "Did the gate land?" reads the gate's FIXTURE ENUMERATION, not its
  assertion. → M-1.
- **L-023** "N spellings of one concept": read the CALLEE's `return` first —
  N callers re-deriving a fact the callee computed and dropped moves the fix
  one hop up. Grep the DEFAULTS of any boolean claim field (`converged: bool
  = True` is the hardcoded lie's mechanism; the type-level fix reds the tests
  that pin the default — rewrite them, do not delete). The brief's named
  EXEMPLAR is a claim (the "deliberate `=50`" helper defaulted to 4000). One
  grep of `filterwarnings` decides raise-vs-warn. → M-1.
- **L-024** A solver's nesting is a per-ENTRY-POINT fact (two public entries,
  same math, structurally different trees; CP has a per-group inner SN does
  not; MoC's inner is a fixed count). Get the tree by RUNNING: a 30-line
  `sys.setprofile` probe with a WATCH set prints nesting, per-level counts
  and truncation (30 of 30 inners at the cap — invisible to any read).
  Protocol-dispatched calls show 0 `callers`. A status derived from a LENGTH
  is guessing on the boundary (49 both for "exhausted" and "converged on the
  last check"). → M-4, M-5.
- **L-025** "What is in scope at P / could P compute Q": patch the callee,
  read `sys._getframe(1).f_locals`, put the ATTEMPTED computation inside
  with a `try/except` that prints the exception type, run on the
  discriminating configurations (slab / sphere-carrying / 2-D-windowed / LD:
  four answers), `inspect.signature` to bound a hole's reach. A `[M]` on a
  negative claim certifies that SOME measurement answered SOME question —
  re-measure it against yours; the mechanism is question-drift, not time
  (the `articulation` brief line carries it since 2026-09-21). "Can A be reused for B?": hunt the input where they DIVERGE most
  (the shipped Krylov inner put them 10⁶ apart), and read the rhs's
  provenance to the loop that produced it (lagged by one outer). → M-1, M-5.
- **L-026** "Is this degraded condition DELIBERATE?": run the test at the
  healthy setting and re-evaluate its OWN assertion — a non-degeneracy floor
  FAILED at the converging budget, a "pinned" gate survived with 20–170×
  headroom; neither is readable. Ask per SOLVE, not per test (the per-test
  CALL COUNT column predicts which rows hide a second bad solve). Never cost
  a budget raise by the ratio (a 3.5× knob cost 1.6× wall). Harvest the stale
  docstrings the counterfactual exposes. → M-5.
- **L-027** Reconciling a survey: `git diff --stat <baseline>..HEAD -- <the
  PRIMITIVE's module>` and read that diff — a new `__post_init__` guard is a
  new obligation on every producer, invisible to claim-by-claim re-checks.
  A handed count that reproduces under no convention was estimated: print
  several conventions side by side, and check the qualifier ("all in
  tests/") as a separate claim (it hid two `examples/` readers). → M-1, M-2,
  M-3.
- **L-028** A UNITS / RANGE / SIGN / ORDER change to a shared producer: table
  the consumers by the GUARD each sits behind and whether it reads the
  changed quantity, then each unguarded consumer's test count. The guarded
  production path is the SAFE one (it raises); the derivations sibling with
  no guard and no test is where the garbage lands. Check the ORDER (ascending
  → descending negates every barycentric) and a ZERO in the new range (a
  normalising division → `inf`). Note when no test asserts the producer's
  VALUES at all. → M-5.
- **L-029** A "build the capability vs refuse the input" fork: price the
  capability's OUTPUT against the incumbent at equal cost (a marginalised LS
  rule loses seven orders to `gauss_legendre`). "Does the primitive exist?"
  is per TIER: read the candidate's MATCHING/precondition code, not its
  docstring's verb (`quotient` refuses anything that moves a node). A guard's
  QUANTITY may be annihilated by the input class's symmetry (the first
  moment is 0 on every `O_h` rule). Sweep every shipped rule through the
  real constructor and read the traceback FRAME. → M-5, M-6.
- **L-030** Before instrumenting a seam, `grep -n "return <Type>(\|return
  _<tail>("` and patch EVERY site: a "the ONE construction site" docstring
  was false (10 → 24 affected tests, silently). Build the census with two
  counters at different depths so "config built, seam never reached" is a
  row you must explain. `[REMEDIED]` the bypass arms have since been folded
  (`[M]` 2026-09-21: 5 `return _package_solution(`, 0 `return Solution(`).
  → M-1, M-3.
- **L-031** A `file:line` deliverable cites by GREPPING the anchor string,
  never from a `sed` window's offsets (18 of 19 off by 3–10, uniformly);
  batch one `grep -n "A\|B\|C"`; spot-check three before shipping — if one
  is off, all are; prefer the anchor string over the number.
- **L-032** Equation-implementer hunt: the `vv-status rationale` naming the
  answer sits on a NEIGHBOURING label (8 of 17: the test's `pytestmark`
  comment, a sibling's sentinel, another page) — grep the LABEL across
  `docs/` and `tests/`, not ±20 lines. A NEGATIVE sentinel ("carries no
  vv-status because …") is evidence a `vv-status` grep misses. Before
  `NOTHING:identity`, ask what object the LHS is and grep for a routine
  returning it (`mass_matrix`, `balance_residual`). `equation_labels=(...)`
  on case dataclasses are existing declarations. A dataclass field resolves
  `py:attribute:`, illegal as an `implements` source — escalate to the class.
  KIND prior: a page of ALGORITHMIC RULES is ~100 % declarable (MC 22 of
  22), an ALGEBRAIC-LAW page ~50 % (operator algebra 21 of 40:
  `{identity, law, canonical-form}` → NONE; `{typing-rule, definition}` →
  look for a declaration site); the CP page's rationale comments answered
  0 of 15, the tests' `pytestmark` comments all 15.
- **L-033** The CODE may already declare the label (`grep -rn ":label:
  <name>" orpheus/`; `grep -rn "(Eq. <name>)" orpheus/derivations/`), and a
  page may say so. The `@pytest.mark.verifies("<label>")` claimant's BODY
  names the SUT. The principled NOTHING's tell is *independent* / "NOT the
  production" (declaring it would make the gate a self-comparison). An
  equation asserting `A ≡ B` can be HALF-retired: existence-check both
  sides. A `w`-generic primitive legitimately implements several labels —
  say so and let the declarer rule on breadth. Tier the answer (arithmetic /
  generic-at-a-constant / factory / symbolic algebra-of-record).
- **L-034** "Doctrine A replaces B" needs THREE inventories: B's
  implementation; A's existing fragments (the flag nobody reads, the
  refusal, the sibling family's battery, the normalisation — these flip
  BUILD to UNIFY/CONSUME); B's own concessions to A. → M-6.
- **L-034b** (formerly the second "L-034") Plan-vs-HEAD includes
  plan-vs-OTHER-PLANS: `ls .claude/plans/ | grep -i <topic>`, and `git
  ls-files --error-unmatch` on the plan itself (a tracked sibling outranks
  an untracked report). A "promote/reify X" row is a claim X exists as code:
  `git log -S "X"` distinguishes retired from never-existed. Run the
  campaign's own gate suite (`-rx`): the xfail rows are the premise oracle.
  → M-1.
- **L-035** A cold reconciliation opens with `git log --oneline --since=<audit
  date> -- <audited files>` (two docstring-only commits ⟹ hunt what the
  audit never READ: 3 of 6 space files). A claim CONFIRMED at its cited site
  is unconfirmed until the site's CONSUMERS are counted — both cited
  factories had zero production consumers while inline mints carried a
  different axis order: "CONFIRMED (substance) / REFUTED (site)". The
  deliverable's unit is claim → verdict → command. → M-1, M-2.
- **L-036** A "(v1 audit, retained)" section keeps its ORIGINAL vintage: bound
  `--since` per SECTION, never per document (0 commits since the plan date
  proved the staleness INHERITED — the plan was written against a memory).
  Grep "retained", "as in v1", "audit kept". When a plan prescribes a design
  the tree has replaced by a ratified counter-design, the verdict is
  "conflict to adjudicate", with both dates. → M-1, M-2, M-6.
- **L-037** AST walks are filters too: `Assign` ≠ `AnnAssign`, `FunctionDef`
  ≠ `AsyncFunctionDef`, `Name` ≠ `Attribute` — a zero on a known-populated
  file indicts the node-type predicate (X1's positive control, carried by
  the brief line); when it matters, `exec`/import the module and read the
  object. A method on class A and a FIELD on class B can share a name
  (`streaming_terms`): a member census needs the RECEIVER column. → M-3.
- **L-038** A "who consumes X" census runs the AST pass FIRST (`Call`/`Name`
  loads in files whose imports bind the name; string annotations counted
  apart) and reports grep lines only as the PROSE column: 3 of 6 public
  `symmetry.py` functions had 0 executable callers behind 4–5 docstring
  `:func:` lines, and Nexus `impact` inherits the flattery (a `:func:` role
  mints a `references` edge). `dead_functions` on those modules: 8 of 8
  false positives (property bodies, catalogue-dict capture, class-body
  install) — confirm by AST (skill: "candidates, not verdicts"). `__all__`
  vs the AST public set is a one-line diff worth printing. The `instrument-doctrine` brief line carries the AST-first clause since
  2026-09-21. → M-3, M-4.
- **L-039** Nexus `callers` on a METHOD returns `nodes: []` plus an
  `unresolved` block; the `unresolved.count` IS the census (it matched an
  AST `Attribute.attr` count to within one) — pair them, never read `[]` as
  uncalled (skill caveat: empty ≠ dead). An entry-result consumer question
  (`.angular_flux` off `solve_sn(`) needs a two-level resolver (direct
  Assign catches 7 of 62; helper-return and fixture-param the rest). When
  the population that must witness a fix is EMPTY, the witness is a
  deliverable to name. → M-3, M-4.
- **L-040** A retirement's PROSE census sorts each hit into (a) an xref that
  DIES with the symbol; (b) a PROMISSORY present-tense claim ("until the S3
  flip", "stays legal until") — MUST-FIX whichever way the ruling goes; (c) a
  CONVENTION statement, false only if a specific ruling lands. Merging (b)
  and (c) over-prices the sweep 4×. `dead_references` measures the PRESENT
  tree (0 before a retirement is the expected reading). A guard-message pin
  census maps each `match=` to the OPERATOR it exercises. A step CODE is a
  homonym across campaigns ("S3" ×3). A chartered assertion with 0 spellings
  flips "fix the prose" to "the step OWES the gate". → M-3.
- **L-041** An identity predicate spelled with `is` on an OPTIONAL
  constituent is vacuously TRUE on the `None` arm (d=3 meshes with different
  edges compared equal) and FALSE on a same-data rebuild (a fresh adapter per
  call): run it on both, report both as `[M]` with the fixture. `[REMEDIED
  2026-09-12]` `is_same_phase_space` → `same_phase_space` (#459 CLOSED).
  → M-5.
- **L-042** Dataclass introspection lies: a `frozen=True, eq=True` record with
  an array field has a generated `__hash__` that RAISES and an `__eq__` that
  raises on ≥ 2 elements. Only `hash(a)`, `a == b`, `a == a` at runtime
  answer, on an `ng ≥ 2` fixture (`ng=1` is a false green). Report the
  introspected flags only as "what introspection says". → M-5.
- **L-043** A kwarg's move census is an AST keyword census split BY RECEIVER
  (148 sites vs 183 grep lines; 54 `La13511Case(scattering_order=…)` were a
  record-FIELD homonym). Trace the VALUE per entry to where it is last
  transformed (clamped on `SNSolver.__init__`, padded on the adjoint entries,
  raw to DSA): a knob "lives" where it is last transformed, and per-entry
  disagreement makes the move a twin-path repair. Grep the string class
  (`"X"` in dicts, `cfg.get`) separately. → M-3.
- **L-044** A "bit-identical" verdict: state each wall's PREDICATE class
  (band vs `nulp`/`array_equal`) beside it — every eigen fixture is a
  `SAFETY × conv_tol` band; the only `nulp=1` walls cover `L + C`. ULP-probe
  the two spellings at the operator level on a random operand (2–5 ULP
  apart) and name the pre-carve `nulp` anchor that must be captured. Spy
  pitfalls: patch the CALLER module's binding (`solver.py` holds its own
  import; the definer's binding read 0), key the spy on BUILDS not lookups
  (207 vs 1), assert the exact count. → M-4, M-5.
- **L-045** A brief's `Class.attr (file:line)` names a LINE, not an owner: run
  it through the file's class boundaries (`awk '/^class |^    def /'`) —
  `:441/:515` sat in `CollisionCache`, not the cited geometry class. A
  `\.method\(` regex counts docstring prose and misses the duck-typed
  `def method` surrogate (`_SpyStrategy.sweep`) and the `setattr` /
  `patch.object` string spellings: three populations, reported separately. A
  one-off disagreement with the brief is almost always a prose line or a
  surrogate. → M-1, M-3.
