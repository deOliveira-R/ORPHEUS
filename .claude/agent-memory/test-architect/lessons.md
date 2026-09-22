# Test Architect — Lessons (hot digest)

Read whole at the START of every dispatch. ONE imperative per entry, stripped to
the correction. Every measured number, war story and `file:line` lives in
`lessons_archive.md` (sections L1–L86, war stories; L87–L88, cross-campaign
lookup tables for families 6 and 8) — open ONLY the section a pointer names.

**Cite, never restate.** The failure-mode taxonomy (AI modes 1–6, test-design
modes 7–12), the three pillars and anti-patterns #1–#36 are `vv-principles`; the
four evidence statements are `instrument-doctrine`; the retirement census is
`retirement-audit`; the denominator clauses are `plan-authoring` §2; the `-O`
scope and the guard-is-debt rule are `coding-standards`; the canonical
invocation, the tolerance contract and the math-type law rule are `vv-testing`;
the reference inventory and the XS mixtures are `AGENT.md`; the per-carve
recipes are `MEMORY.md` §3. An entry exists HERE only because no rule and no
skill carries it.

**THE SPINE** (standing: `AGENT.md` §0.5 / §0.6 / §1.5). A plan is done not when
the tests pass but when, for EVERY gate: (a) a named mutation reddens it under
`python -O`; (b) the reference is structurally INDEPENDENT of the SUT; (c) the
regime ACTIVATES the term the bug lives in.

**Maintenance.** A new lesson adds its RULE here with a `→ LNN` and its war
story as a NEW archive section — never a campaign-named block appended here. An
entry recording a GAP gains a LANDING note the day the gap closes: `L3`'s "no SN
MMS exercises `q.boundary ≠ 0`" was true when written, the §4.6 family closed
it, and the stale entry then generated a whole phase brief for work already
done (`L40a`).

---

## 0. Eight meta-lessons — most entries below are an instance of one

**M1 — CONSTRUCT THE SUBJECT AND MEASURE THE ASSERTED QUANTITY BEFORE DESIGNING
THE GATE.** Not read the design, not read the docstring: build the object the
carve will produce, print the field the gate will assert, evaluate the dispatch
predicate on a production instance, probe the FACE the consumer actually reads.
Half of all findings are that the gate is unwritable, vacuous, or already true.
→ now in `AGENT.md` §1.5, a standing directive (2026-09-21).
→ `L41b`, `L43h`, `L58f`, `L67a`, `L77a`, `L77f`, `L79b`, `L80f`, `L81c`,
`L83a`, `L83c`, `L84i`, `L86d`

**M2 — READ EVERY NULL THROUGH FOUR HYPOTHESES BEFORE "THE GATE IS BLIND":** the
mutation was INSUFFICIENT · a Pattern-2 TWIN predicate survives and guards it ·
the geometry or type ANNIHILATES the degree of freedom · the INSTRUMENT never
installed. The last is cheapest to check and the most common.
→ `L49c`, `L68a`, `L68c`, `L73j`, `L75a`, `L79f`, `L81d`

**M3 — COMPUTE THE STABILISER OF THE WHOLE GATE SET AT DESIGN TIME** (`vv` Mode
12), and remember there are TWO sides. OPERATOR side: `[G, Aᵀ] = 0`. SPACE side
(the dual nobody had written): a rank-1 point axis makes every weight a
one-element array, i.e. a SCALAR, and a scalar `G` commutes with everything. A
ratio annihilates uniform scale; a spectrum, similarity; a palindromic rule, its
own reversal; a symmetric generating rule, half the permutation group.
→ now `vv-principles` Mode 12's space-side check (2026-09-21).
→ `L43e`, `L47a`, `L47e`, `L58b`, `L59a`, `L61a`, `L84i`

**M4 — EVERY CLAIM ARRIVING IN A BRIEF, PLAN, CHARTER, DOCSTRING OR DESIGN MEMO
IS A HYPOTHESIS TO RUN.** Refuting the optimistic premise with a MEASUREMENT,
before the ink dries, is the highest-value output of a proactive dispatch; state
the refutation IN the plan so the implementer ships the achievable carve.
Enumerate the population yourself (`dir(module)`, `__subclasses__`,
`dataclasses.fields`), never the brief's list.
→ `L10`, `L40a`, `L43a`, `L63d`, `L64g`, `L74b`, `L80b`, `L81g`, `L86a`

**M5 — THE TREE MOVES UNDER YOU.** The checkpoint protocol (`git status
--porcelain` and `git log` at the START and the END; a pristine copy before the
first mutation, `diff -q` after) is `process-discipline` ("Trust git",
"Mutation-testing an uncommitted file"), cited, not restated. What this tree
adds: `ls`/`shasum` the target before any pre-carve claim and bracket every
measurement with the mutated thing's state, because a concurrent write can
DEMOTE your measurement to a value-compared-with-itself; another agent's
deliberate mutation is indistinguishable from a production bug (re-read with
`inspect.getsource` before reporting one); a carve can land twice, three times,
or mid-sentence, and then every gate becomes MEASURED rather than predicted.
→ `L28`, `L33`, `L43a`, `L44h`, `L72a`, `L75`, `L84a`

**M6 — A UNIVERSAL OWES A DENOMINATOR MEASURED OVER THE SHIPPED POPULATION, AND
THE AXIS MATTERS MORE THAN THE NUMBER** (X2, `plan-authoring` §2). Axes each
wrong at least once here: per ARM, per MEMBER of a union, per BRANCH of a
dispatch, per CONSUMER, per CALL SITE, per ROW of a parametrize, per FAMILY of
the corpus, per SUB-FAMILY of a tolerance. A corpus uniform in the
discriminating field leaves an arm witness-less; a single *member* can make an
invariant unspellable by construction.
→ now `plan-authoring` §2 AXIS-CHOICE (2026-09-21).
→ `L59b`, `L60e`, `L64f`, `L65d`, `L66d`, `L70e`, `L71d`, `L74b`, `L76f`,
`L77b`, `L86e`

**M7 — WHEN TWO ROUTES ARE NUMERICALLY IDENTICAL, NO VALUE GATE CAN SEE THE
CLAIM.** The instrument is then IDENTITY (`is`), ROUTE (swap the old owner for a
decoy) or COUNT (a call spy). Reach for one the moment a carve's own
justification is "bit-identical by construction".
→ `L64a`, `L65c`, `L66b`, `L77e`, `L77h`, `L82c`, `L83f`, `L85a`, `L86c`

**M8 — A CARVE CHANGES WHAT THE SURVIVING GATES MEASURE**, four ways, all
silent: DEMOTE (two sides become one object), PROMOTE (the gate got stronger and
its docstring still advertises the weak claim), DIE (it can no longer CONSTRUCT
its subject — delete it, never repair it by passing the new argument), INVERT
(it now pins the degradation as the contract). Inventory the survivors BEFORE
the carve and re-pose them in the SAME commit.
→ `L34e`, `L36`, `L37`, `L61c`, `L66f`, `L75e`, `L83d`

---

## 1. Gates that cannot red — the shapes `vv` Mode 8 / anti-#17 do not carry

- **⛔⛔ RETIRING A DEFERRAL: never read the reason string — RUN the row.** The
  fix an exclusion CITES is usually not the fix that healed it. Method: a `-p`
  plugin no-op'ing `pytest.xfail` that ASSERTS its own installation, one run per
  row; attribute the cure afterwards, and say UNDISCRIMINATED rather than pick
  between two candidates. Riders: a reason's named BUDGET may not be a live
  knob; an UNCONDITIONAL stub whose body is only the `pytest.xfail` call must be
  SUPPLIED a body with exactly one failable statement; the healed row's doc
  claims go present-tense-FALSE. Sharpens `vv` Mode 8(9). → `L45`
- **⛔ The xfail family's four silent failures.** (a) A flip-edit must touch a
  statement whose VALUE the production change determines — diff the xfail body
  against its own flip-proof; textually equal ⟹ ceremony `L33`. (b) A marker
  SPLIT (`@xfail` → `pytest.param(marks=)`) loses `strict` silently, and
  `--collect-only` and `-rx` are BOTH blind (`[M]` 2026-09-21: `pyproject.toml`
  still has no `xfail_strict`); the 5-line permanent catcher introspects
  `pytest.param(...).marks` at import `L61g`. (c) A strict xfail flips only on
  XPASS, so an API that lands WRONG leaves it `xfail` and the suite green — pair
  it with a RECORD row (green today, designed to RED at the carve, DELETED not
  repaired) `L82e`. (d) An IMPLICATION row passes VACUOUSLY when nothing
  satisfies its antecedent, which a pre-carve tree guarantees: owe it a
  non-vacuity guard and a strictness leg `L82d`.
- **⛔ Refusal gates: the message IS the gate.** `pytest.raises` without `match=`
  legs KEYED to the argument that triggers it is teeth-less, and a blanket "the
  message names both completions" pins a FALSE reason on the row whose defect is
  different. A `-> T | None` refusal verb collapses N guards into one value:
  isolate from the INPUT side (an input passing all guards but one) and mint
  DISJOINT fragments, asserting the disjointness once. When the fix is "raise
  the project's TYPED error instead of the builtin", `pytest.raises(<builtin>)`
  is green BEFORE and after (`BoundaryError(ValueError)`) — name the SUBCLASS
  and require the pre-fix state to be RED. → `L31`, `L35f`, `L36f`, `L75e`
- **⛔ A new guard wired AFTER an existing one: the earlier guard's inputs are a
  DISCRIMINATION row, not a negative row** — assert, on that same input, old
  fragment PRESENT + new fragment ABSENT, which pins the WIRING ORDER. When a
  new law and an existing clause share one `__post_init__`, the old clause's
  historical witness violates BOTH, so each owes a DISCRIMINATING input and the
  ORDER owes its own row. → `L43c`, `L73g`
- **⛔ A guard's REACHABILITY is decided by whatever runs before it.** Enumerate
  the guards ABOVE a refusal and the PRODUCERS of its predicate before
  re-wording it; a refusal no shipped input reaches is RETIRED, not re-worded.
  Put a typed guard at the TOP of the promoting classmethod — the caller's own
  attribute reads raise a frame earlier, so the error TYPE is the placement
  gate. The audit's fifth search: a retirement can ORPHAN a guard by removing
  its only public ROUTE while guard and witness both survive. → `L60f`, `L79d`,
  `L81c`
- **⛔ Retiring or shedding a FIELD disarms every guard that KEYS on it** — grep
  removed names as GUARD PREDICATES (`is None`, `getattr(…, default)`,
  `hasattr(`), not only as reads, and land the re-key with its witness in the
  same commit (`vv` #28's temporal twin). Sibling: an attribute→property
  conversion kills every `hasattr(Class, …)` PREMISE, and a READ census cannot
  find a premise. → `L60c`, `L62b`, `L63c`
- **⛔ Two quantifier traps.** A ∀ over a per-element predicate is UNGATE-ABLE
  when no factory produces a MIXED input; the fix is architectural — return the
  offending POSITIONS (`vv` anti-#14) `L43d`. And count the rows that REACH the
  assertion, not the rows that exist: a guard-clause early return is anti-#20
  wearing a guard clause `L43f`.
- **⛔⛔ `warnings.warn(stacklevel=N)` is a claim about EVERY call site's DEPTH
  and NO message gate can see it**; worse, the obvious gate is blind to half the
  class ("the attributed file is outside the package" reds for `→2`, stays green
  for `→4`). Ship TWO legs — `not is_relative_to(pkg_root)` over every entry,
  PLUS `w.lineno == inspect.currentframe().f_lineno + 1` above a DIRECT call —
  and never gate the literal `stacklevel == 3`. The structural fix beats the
  gate: hoist every emission into the public entry. → `L46b`, `L46c`
- **⛔ A TYPE-ANNOTATION widen has NO runtime witness** ("hand it the wider type,
  assert it constructs" is green before and after). Gate it with
  `tests/test_pyright_ratchet.py`, and SAY which done-when items are grep
  OBLIGATIONS rather than gates, or they read as covered. → `L59d`
- **⛔⛔ A flagship NUMERICAL gate can be a THEOREM with no reachable falsifier —
  grep a charter for *"for the wrong reason"* / *"the same result today"*.**
  Repair: gate the theorem's PREMISE, which IS red-capable; keep ONE corollary
  row labelled claim-kind THEOREM carrying the blindness table; name the
  pre-existing `vv` #19 control as the only loaded partner; do NOT manufacture a
  wrong-structure control the production type refuses to construct. A
  reverse-composite law (`(RKM)† = M†K†R†`) is a theorem of the metric adjoint
  and cannot gate the faces it is built from — an identity holding for ARBITRARY
  factors is `vv` anti-#24(d) whatever it is named. → `L61a`, `L67a`, `L73a`
- **⛔ A gate that builds its reference THROUGH the object under test sees only
  self-consistency** — a null vector from `svd(A)` makes "blind to `ker A`" a
  fact about the FACTORISATION. Re-pose onto the MEASURAND production reports.
  → `L35k`, `L49a`
- **⛔ A pin naming a "legacy"/"reference"/"adapter" counterpart: two probes, in
  order, before any battery.** (1) Is the other side literally the SAME OBJECT
  (`is`, five seconds — a shared `cached_property` made one leg `array_equal(x,
  x)`)? (2) Garbage the ONE shared producer in EVERY module binding. Then
  RE-SCOPE, never delete. (`retirement-audit` D.14 at the fixture tier.)
  Companion: a docstring naming "the surviving pins" is a CLAIM to measure,
  including one you wrote an hour ago — and an L0 identity that "covers" a term
  may RECOMPUTE the production array instead of reading it, pinning the LAW
  while blind to the ARRAY. → `L34e`
- **⛔ Retiring a runtime guard that had NO negative test makes its replacement's
  teeth NET-NEW, not migrated** — grep `pytest.raises(match=<guard msg>)` before
  crediting a mechanism-swap as behaviour-identical, and write the negative test
  the guard never had. Production guards in a step's blast radius routinely have
  ZERO witnesses; grep each `raise`'s shortest distinctive fragment. → `L4`,
  `L58g`, `L83j`
- **⛔ When a defect was closed STRUCTURALLY the obvious mutation reds NOTHING**
  and the `catches` marker looks unearned: the TYPE refuses the bad value one
  frame in. Target the type's invariant, say so in the marker's docstring (a
  stronger claim). General form: a law landing as PREVENTION-BY-CONSTRUCTION
  reds 0 and makes builder-level mutants UNINSTALLABLE — only a DIRECT
  construction supplies the witness, and an uninstallable arm IS a finding whose
  content is the refusing guard's NAME. → `L41d`, `L72d`, `L73k`
- **⛔ Mutations that red by RAISING attribute nothing** — a refusal making an
  old spelling UNCONSTRUCTIBLE turns the "revert" arm into a crash arm; an
  out-of-range permutation is out of range BY A THEOREM. Always ship the
  attributable twin: an IN-RANGE, in-class mutation (right set, wrong
  assignment) that reds by COMPARING. A mutation INAPPLICABLE to the fixture's
  shape is not "no teeth". → `L25`, `L31`, `L41c`, `L42a`, `L70d`
- **⛔ Constructing a break-exactly-ONE-invariant mutant is a design problem** —
  `np.roll(arange(N),1)` breaks measure AND sign AND involution, and no ODD
  cycle can isolate an involution. Carry the fixture PER ROW. → `L31`
- **⭐⭐ THE ROUTE GATE, for "consumer X now reads from OWNER A, not B" — no
  value or identity gate states it.** Pose, SWAP the old owner's object for a
  decoy, require the answer UNMOVED (or MOVED, for a re-point). Three traps:
  mutating ONE consumed surface certifies one route; a DRIVER that re-poses
  internally measures nothing; a surviving cache MASKS the swap, so the gate
  needs an ACTIVATION leg. The DECOY must clear the PRODUCTION ADMISSION GUARDS
  of the arm the gate lives on, not merely discriminate — print its
  discriminating array first. → `L64a`, `L65c`, `L66b`
- **⭐⭐ When two production routes are NUMERICALLY IDENTICAL the instrument is a
  CALL COUNTER**, and its first red is a 0-vs-N contrast on the SAME suite:
  install at `pytest_configure`, count the SPECIFIC verb, assert the SIBLING
  counts unchanged so the row is attributable. A route claim is about the CALL,
  not the callee — assert the MECHANISM. And a spy on a SHARED verb observes the
  CALL ARGUMENT, not the value the verb serves: the right instrument only when
  the claim is about the CALLER, and that sentence belongs in the docstring.
  → `L77h`, `L82c`, `L83f`
- **⭐⭐ A shape-keyed SCANNER used as a ruled row's predicate is a FILTER.** Its
  flip-proof must plant a DECLARED member (an `object.__setattr__` staple is
  invisible to a `dataclasses.fields` walk), and it ships with a planted-member
  POSITIVE control, a NEGATIVE control on a member-less neighbouring shape (else
  it over-matches), and both scanners HOISTED to module scope. → `L86b`
- **⛔ A non-vacuity leg must be EVIDENCE, not an assertion inside an xfail** (an
  xfail hides ANY failure, `vv` Mode 8(4)): wrap the call in `try/except`, fold
  the exception into the message, leave exactly ONE failable statement. Sibling:
  a pass-through row asserting EQUALITY is usually unreddenable — `assert out is
  q` gives it teeth, and the tell is that no arm in your OWN battery touches it.
  → `L78p`, `L86c`
- **⭐ Ship the arms designed to go GREEN, and read the ones that do not.** A
  DECLARED PARTIAL NULL arm is the only instrument that can state a flagship
  gate's own Mode-12 blindness. A DECLARED-BLIND arm that REDDENS is a finding
  about the gate SET — read its red SET, which can PARTITION the suite as no
  green arm could. Rows NO arm reds may be DECLARED reference-side controls: say
  so, because their green is the LICENCE to read the neighbours as coverage, not
  coverage. → `L72f`, `L73i`, `L75a`, `L78f`, `L81g`, `L83h`
- **⛔ A declared blindness must name the RIGHT symbol and stay inside the
  problem's CONVERGENT regime**, or the fixture's convergence guard fires first
  and every red reads *"did not fully converge"*. Two more: a defect with TWO
  ENDS needs TWO arms, and the arm you write covers the end you were looking at;
  and a solver entry's STRATEGY DEFAULT decides which production branch a whole
  module ever poses — `inspect.signature` the entry and enumerate the defaults
  before claiming branch coverage. → `L78l`, `L78m`, `L78o`
- **⛔ When a "remove the step" mutation is NULL, try "corrupt the step"** —
  idempotence at a fixed point makes REMOVAL invisible and CORRUPTION loud. A
  step measured inert is a NUMBER in the commit message, never a gate. And a
  "complementary pair" of guards can be complementary in ONE variable and both
  silent in another: draw both complements and check they COVER. → `L78b`,
  `L78e`
- **⛔ Two defect classes with the IDENTICAL law residual need a STRUCTURAL leg**
  — ship the partition as two separately-messaged assertions. And a `(1+ε)`
  VALUE mutation cannot red an `is` row; a battery predicting it predicts the
  impossible. → `L84j`, `L84l`
- **⛔ The identity family: three legs that are false reds waiting.** (a) A
  done-when spelled `is` on a VALUE type — `__post_init__` normalises, so the
  constructor returns a FRESH instance; gate a value-merge with
  `==`/`hash`/`name`/`repr`/container-dedup and say `is` is NOT asserted. (b)
  `hash(a) != hash(b)` is NOT a legal "these differ" leg (a frozen dataclass
  hashes the FIELD TUPLE, not the class) — assert separation through the
  CONTAINER, and run a constant-`__hash__` arm before ANY identity carve: its
  red set is the re-pose list. (c) A value-MERGE collapses roster denominators
  silently — assert `len(set(roster)) == len(roster) == N` and grep the roster's
  docstring for the sentence the merge falsified. → `L70a`, `L72b`, `L72g`,
  `L79e`
- **⛔ Price every guard at EACH TIER, and ask which arm produces the headline.**
  A guard load-bearing at its API tier can be INERT where its end-to-end test
  lives (a later stage refuses first) — write the inert tier into the docstring.
  A design's HEADLINE consequence can be its least falsifiable one: ask whether
  its arm READS THE DATA. And before writing a control that separates two
  groups, ask whether their INVARIANT RINGS differ — coinciding rings make the
  control unwritable at every fixture, and the honest deliverable is the
  measured inertness as a NAMED blindness row. → `L70b`, `L71c`, `L71f`
- **⛔ Before any linearity / homogeneity / additivity row, measure `|Op(x)|` on
  a random `x` and require `> 0` as a committed ACTIVATION leg** — where the
  operator collapsed to the zero morphism both sides are structurally zero and
  no input can red it; the honest gate is then the STRUCTURAL claim on that
  fixture plus the linearity row on one where the operator is non-trivial.
  → `L40c`
- **⛔ REPAIRING a decayed gate is a different design problem from writing one.**
  Re-pose onto a REGIME-INDEPENDENT mechanism, never drive the fixture back into
  the regime; check reachability BEFORE trying to reach it; never compute the
  reference with the routine that ESTABLISHES it. The acceptance measurement is
  the PER-GUARD red table with every residual miss categorically out of scope
  and said so; re-run the AUDITOR's own harness, never a re-implementation.
  → `L28`
- **⛔ A design memo's HAZARD PROSE is a claim — run it.** A characterization
  test freezing a no-guard ruling asserts CONSTRUCTIBILITY only (one positive
  leg, no negative) and says why; the ruling's justification sentence must not
  reach the constructor docstring, or it reads as licence to add the forbidden
  guard. → `L64g`
- **⛔ A gate can stay green while its REASON becomes false** — when a phase
  falsifies a structural claim, grep the claim's WORDS in `tests/`, not its
  symbols (`retirement-audit` A.7), re-scope in the SAME change, and give the
  new structure its own positive gate. → `L33`

- **⛔ A DIMENSION count (SVD nullity, rank) is blind to a kernel ROTATION.**
  When `A = A_RR ⊕ 0_K` with `A_RR` nonsingular, a defect that fills `A[R,K]`
  leaves the nullity EXACTLY unchanged (the kernel becomes `e_t − A_RR⁻¹A_RK e_t`),
  so a law stated as `dim ker A = … + |K|` stays green. Gate the OBJECT: `A e_t == 0`
  and `(Ax)[K] == 0` bitwise, one leg per side (forward rows, transpose rows), each
  with its own tooth. And a reciprocity gate on a singular metric sees only
  `A[range, ker]`: writes into kernel rows and kernel→kernel maps are its
  stabiliser. → `L90`

## 2. Harness discipline — the instrument lies before the code does

`vv` anti-#17's nine checks and `instrument-doctrine` X1 are the rule. Below:
only the ORPHEUS mechanisms they do not name.

- **⛔⛔ NEVER quote a NEXUS-derived per-node test count as coverage — in THIS
  codebase it measures the RESOLVER, not the suite.** Static `callers` missed
  **217 of 229 (94.8 %)** of the tests that execute `OperatorSum`; **21.3 %** of
  `calls` edges tree-wide point into `unresolved`. ⭐ The severity is a
  CONSEQUENCE of Cardinal Rule 2 — `coding-elegance` Pattern 1 spells every
  operation as a dunder on a domain type, exactly what the resolver cannot
  follow, so **the better the architecture gets the blinder the call graph
  becomes.** The repair is the runtime overlay. The only tell that caught it was
  IMPLAUSIBILITY, not the instrument. → `L55a`, `L55b`
- **⛔ Make the harness ASSERT its own installation; a banner nobody reads is not
  a check.** Every plugin `raise`s unless it rebinds N of N symbols, and the
  banner COUNT is grepped into the result line. Five measured false-"0 caught"
  verdicts, all flattering: a shell loop that DROPPED the `-p` flag; a loop that
  TIMED OUT (budget from the MUTATED cost — garbage destroys convergence); an
  arm naming the wrong MODULE for a method defined on a BASE (`rc=3 / FAILED=0 /
  banner=0`; against the base, 41 reds); a wrapper defeating
  `inspect.signature`, reddening a helper for INSTRUMENT reasons; a census
  plugin rebinding only the DEFINING module while tests import the package
  re-export (rebind every `sys.modules` entry whose attribute `is` the
  original). → `L46e`, `L79f`
- **⛔ Give every arm a BITE CHECK comparing the mutant's ANSWER to the honest
  one at the load-bearing input, and evaluate the honest value BEFORE the
  patch.** Four ways a bite check lies: a per-instance MEMO warmed by an earlier
  solve masks the mutation and reads a plausible bit-identical GREEN (install at
  `pytest_configure`, before any object exists); a captured "honest" callable
  that dereferences the patched name by LATE BINDING compares mutant with
  mutant; `assert SUT is mutant` proves the REBIND, not the bite; and `apply =
  _apply_impl` is an ALIAS, so rebinding `_apply_impl` changes nothing the alias
  sees — wrap `cls.__dict__[verb]` and call `descr.__get__(self, cls)(x)` so
  `singledispatchmethod` still dispatches. → `L68b`, `L73j`, `L74i`, `L77e`
- **⛔ A "deleting X reds 0 of N" measurement is VOID when X is imported at
  MODULE SCOPE on a conftest's import chain** — one notch past `vv` Mode 8(3):
  `rc=4, 0 collected, 0 ^FAILED AND 0 ^ERROR`, and
  `--continue-on-collection-errors` does not help. One `grep -rn "import X"
  orpheus/` answers it first; the honest instrument is the IN-CLASS rebind (`vv`
  #18). A class DELETION's collection-killers are MORE than an audit names: a
  module-scope attribute read AND a construction inside a module-level
  `parametrize` list. → `L67i`, `L68a`
- **⛔ MY OWN parametrize list has called production at MODULE SCOPE three times**
  (`vv` anti-#17(c)), twice recorded and committed again. Parametrize over a
  LITERAL label tuple, build in the BODY, and gate the literal against the
  producer inside a body — which also gives the population a falsifiable
  denominator row. → `L47d`, `L81e`
- **⛔ In a mutation plugin a PRECONDITION is the arm's FIRST statement and
  raises a distinct `Uninstallable`, never the bite's `RuntimeError`** — a
  partial install under a failed precondition is worse than a crash. Check every
  arm against the CURRENT tree: a prior carve can dissolve a mutation's own
  distinction. → `L74i`
- **⛔ Build a source-mutant by TRANSFORMING `inspect.getsource`, never by
  hand-copying** — a hand copy is a twin path that drifts, and a `str.replace`
  whose target is ABSENT can `raise`, which makes the instrument assert its own
  installation. Smoke-test each mutant's OUTPUT before the battery. ⚠
  `textwrap.dedent` strips FOUR spaces from a method's source, so a target
  copied at class indentation never matches; and a re-typed mutant's re-worded
  `raise` reds a `match=` gate for a reason the mutation is not about.
  → `L44i`, `L75d`
- **⛔ IMPORT-CHECK the campaign's own mutation harness at the LANDING tree
  before planning a battery on it** — a campaign that retires symbols breaks its
  instruments by MODULE-SCOPE BINDING, silently. One `hasattr` is the whole
  check; repairing the harness is part of the commit, since no negative verdict
  is trustworthy until its control passes. → `L44d`, `L84a`
- **⭐ Prefer a monkeypatch-only / `dataclasses.replace` battery: crash-safe BY
  CONSTRUCTION**, strictly stronger than copy-aside + `diff -q`
  (`process-discipline` § "Mutation-testing an uncommitted file"). Mutate a
  frozen VALUE by rebuilding it, never production. Measure "before" with `git
  show HEAD:<file> > <tmp sibling>` so the sibling runs under the SAME mutation,
  collection and fixtures. → `L28`, `L63h`, `L84j`
- **⛔ MEASURE the reachable subset before rationing the battery** — an
  over-stated cost silently shrinks it, the same loss as a blind gate and harder
  to see (one directory: 329.66 s whole vs 2.72 s for the four files a carve
  could reach). Build the attribution scope FROM the positive control's red set;
  re-derive the scope from THIS phase's cone (an inherited battery is scoped to
  the PREVIOUS phase's blast radius); and let the EXCLUDED numbers justify
  themselves in the plan, or an excluded directory reads as an oversight.
  → `L36e`, `L60h`, `L62a`, `L62f`, `L77i`, `L86g`
- **⛔ Compare a superset battery PER ARM, never as a union** — two arms can
  PARTITION the corpus by geometry, so a union hides one side's regression.
  Re-running a prior phase's battery pins a DIFFERENT claim here, the INVARIANCE
  witness: the red SET must be EQUAL per arm, since a same-sized-but-disjoint
  set means the re-source landed on a different instance. Report the split NEW
  vs PRE-EXISTING per arm, never a total — the arms with ZERO pre-existing
  catchers are the headline. → `L64e`, `L66c`, `L66j`, `L70e`
- **⛔ A red set entirely INSIDE the new gate class is not automatically `vv`
  #17(e)'s "mirror, not a gate"** — discriminate by asking whether the consumers
  EXIST and are blind FOR A STATED REASON, told apart by the PRE-CARVE consumer
  census, never by the count. → `L73l`
- **⛔ Attribute an out-of-scope red by AUDITING THE DIFF for arithmetic**, not
  by re-running before/after (which this tree forbids): `git diff -U0 orpheus/ |
  <strip comments> | grep -E '=|def |raise |return '` is the complete
  added/removed CODE line list. And a pre-existing red must be CHARACTERISED,
  not counted — "1 pre-existing red" yields nothing; "the GL case, 1 ULP, on the
  one rule the change does not touch" yields a `-k` filter and a do-not-absorb
  instruction. → `L42g`, `L44l`
- **⛔ Run PYRIGHT over your own new TEST module** — it catches the elegance
  defect, not just the type (a string-tag parametrize with an `if entry == …`
  chain and a `**kw` splat gave 24 errors, where the `# type: ignore` reflex
  would have hidden a real `coding-elegance` anti-#4). Measure the COMMITTED
  file too: a mis-placed `# type: ignore` hides there. → `L44k`
- **⛔ A numeric table in MY OWN plan is an `[M]` claim — the obvious
  continuation of an integer sequence is NOT a measurement.** Compute the
  extension in the same probe that produced the committed rows, or mark the row
  a placeholder. → `L42d`
- **⛔ Census discipline, five ORPHEUS shapes.** A §6b census's own POSITIVE
  CONTROL is what catches it, and the wrong answer points the FLATTERING way
  (`^\s*Quotient\s*\(` returned 1; `\bQuotient\s*\(` returned 10). Two
  overlapping predicates in a brief are NOT two work items — compute the UNION
  and print both set differences. Triage a CONCEPT grep by MEANING (one word
  named three unrelated things, one load-bearing). Split a retype census by ast
  CONTEXT (Load / Store / keyword / Subscript) — an inherited "38 reads"
  contained 39 loads and NONE of the 18 STORES. A class-NAME census misses a
  SUBCLASS inheriting `__init__` and a `Base.create(**kwargs)` registry call:
  resolve bases transitively. → `L58e`, `L66c`, `L73d`, `L73e`, `L76i`
- **⚠ Two shell mechanisms that read as "the battery found nothing" rather than
  "the battery ran nothing".** zsh does NOT word-split an unquoted parameter
  expansion, so `pytest $SCOPE` passes ONE argument and collects 0 — use
  `${=SCOPE}`, an array, or a driver taking `"$@"`, and print the collected
  count. And a `nohup … &` chained AFTER an `until` loop in one tool call never
  launches when the call is killed at its timeout, while `pgrep -f <script>`
  then matches the dead shell's own heredoc and prints RUNNING — launch long
  jobs with `run_in_background: true` and ONE command. → `L63h`, `L68f`, `L80h`
- **⛔ A two-sided JOIN needs a denominator assertion on BOTH sides; the
  unasserted side is the one that fails** (a clean, confident `JOIN RATE = 0.0 %`
  because one side was JSON-encoded). Validate the DECODER as well as the
  filter, by cross-checking against an independently-vocabularied count.
  → `L46e`, `L55e`
- **⭐⭐ A §6b table is a MEASUREMENT, not a reading: wrap the method in a `-p`
  plugin, return the honest answer, and record `(test id, support, before,
  after-shadow)` over a real suite run** — 1636 calls / 74 tests → exactly ONE
  verdict moves, where a grep returns 61 sites and cannot say which. When a
  concurrent carve holds `orpheus/`, SNAPSHOT it (`git archive <HEAD> orpheus |
  tar -x`) and SHADOW the design outside the package (⚠ the editable install's
  MetaPathFinder beats `PYTHONPATH` — strip it from `sys.meta_path`); the
  shadow's validity control is reproducing the shipped answers EXACTLY on every
  input the design does not touch. When the API does not exist yet, the runnable
  dry-run is a SHIM. → `L71a`, `L71b`, `L73c`
- **⛔⛔ Adjudicate a proposed CONSTRUCTION GUARD by INSTALLING it as a plugin
  and counting reds — a per-INSTANCE census (`vv` #29), never a static site
  count.** One charter's guard destroyed 250 of 845 rows (unrunnable, not weak);
  the alternative was live on 18.8 % of constructions and raised 0 times, so it
  had no witness anywhere. A SITE census counts call LINES and understated the
  inertness 4×; put the INSTANCE-tier number in the guard's docstring. → `L61b`
- **⛔ Run the mutation BEFORE writing rewire prescriptions** — a per-gate
  claim-class verdict guessed from reading is wrong in both directions (25 of 33
  reds were a family the brief never listed; a LISTED gate red 0 because it
  feeds the helper's output to the SUT, `vv` anti-#22). Run the teeth harness
  over your OWN new module before delivering it. → `L34d`, `L63d`
- **⚠ After adding a field to a type, grep the tests for REFLECTION walkers**
  (`vars(`, `asdict`, `fields(`) — a walker over arbitrary objects sweeps the
  new field's arrays into an unrelated count and reddens for the wrong reason.
  → `L65h`
- **⚠ `ast.col_offset` counts UTF-8 BYTES; `ast.get_source_segment` is the only
  safe reader — and quantify such a hazard against the CORRECT implementation,
  never a proxy** (the proxy said 128 of 128 spans corrupt; the honest
  instrument read 0 of 128). A latent hazard with 0 witnesses is a RIDER, not a
  defect. → `L74g`
- **⛔ `dead_references` on an UNCOMMITTED working tree reports GRAPH staleness**
  — settle it with a control-validated grep, not by repairing; re-run after the
  next `sphinx-build`. → `L78n`

## 3. Config blindness — the ORPHEUS fixture-fact inventory

Generic rule: `AGENT.md` §0.6, `vv` §H2 / anti-#3 / anti-#4 / Mode 7. Below: the
builders and corpora that SILENTLY null a channel. Each is a measurement with a
shelf life — check it against a concrete row before trusting a green.

- **`make_mixture` nulls TWO channels**: `sig_2` defaults to all-zero and there
  is NO `sig_l` parameter (it hardcodes `SigL = zeros(ng)`), so any (n,2n) or
  (n,α) term is identically nulled and a "balanced" fixture built through it is
  IMBALANCED by exactly `sig_l`. Build `Mixture(...)` DIRECTLY. → `L1`
- **Every `xs_library` mixture (A/B/C/D × 1g/2g/4g) ships `Sig2 = 0` and
  `len(SigS) == 2`.** So any (n,2n) leg there is vacuous (manufacture it —
  `Mixture.Sig2` is a Legendre STACK), and `min(L, len(SigS)−1)` maps every
  request ≥ 1 to 1, so a row spelled "P2 vs P3" asserts `1 != 1`. The only
  non-vacuous abstract-library truncation pair is `(0, ≥1)` (a real 476 pcm
  discrimination); the 421-group library is the other option. Assert the premise
  IN the row. → `L78d`, `L82a`, `L86g`
- **⛔ A manufactured cross section not balanced into `Σ_t` makes the reported φ
  differ from `∫ψ dΩ` by an EXACT GLOBAL SCALE — and the damage is that the L=0
  CONTROL reds too, so the gate attributes nothing.** Use the `tests/cp`,
  `tests/mc` house spelling `sig_t = sig_c + sig_f + rowsum(sig_s) +
  rowsum(sig2)`. → `L78d`
- **Only mixture A is fissile, and EVERY fissile 0-D mixture is SUPERCRITICAL**
  (`k_inf = 1.5 / 1.875 / 1.4878`); B/C/D give `rank(A⁻¹F) = 0`. Hence a
  `solve_sn` eigen snapshot on a moderator mixture is `k = 0/abs → nan`, a
  silent dead test (reformulate as fixed-source against `φ = (diagΣ_t −
  Σ_s0ᵀ)⁻¹Q`), and a "subcritical infinite medium" row is unreachable as shipped
  — manufacture it with `replace(mix, SigP=0.4·SigP, SigF=0.4·SigF)` → `k = 0.75`
  EXACTLY (k is linear in the scale because `F` is rank-1), with 0.6 flipping
  the sign as the refusal control. → `L7`, `L84d`
- **`A − F` is NEVER positive-stable on the SN composite** (8 eigenvalues at
  exactly −1.0, the TRACE rows, at every k) while the BULK block is. A predicate
  spelled "positive-stable" refuses every SN problem; the honest one is
  `ρ(A⁻¹F) = k_eff < 1`. → `L84e`
- **The SN operator fixtures carry `placeholder_materials`** (SigS / χ / νΣf all
  zero) ⟹ `F` is the ZERO operator and its reciprocity row is `0 == 0`. Every
  `tests/sn/architecture/_config` mesh is NON-FISSILE (`solve_sn` raises
  *"leakage scale bridge is degenerate"*) and reads `|N2N·x| = 0`. "Reuse the
  existing fixture" is a hypothesis: measure it, and record WHY each neighbour
  was rejected. → `L26`, `L83`, `L84j`
- **A cylindrical `SNMesh` admits only CARRYING quadratures** (15 of 15
  `folded_product`; 0 of 20 `product`/`lebedev`/`level_symmetric`). So a
  folded-vs-UNFOLDED equivalence gate is UNWRITABLE there, every admissible rule
  IS the σ_y quotient, and the σ_y half of the cylinder's declared symmetry is
  Mode-12 blind BY CONSTRUCTION — the honest gate asserts the structural
  commitment (the solver refuses every unfolded rule). → `L68e`, `L75c`
- **A 1-D rule carries `μ_y = μ_z = 0` on every ordinate** (ERR-080), so any
  σ_y/σ_z evenness leg there reads `0.0` with IDENTITY permutations: owe every
  such leg a `perm moves k/N` VACUITY guard. Read a fold's quadrant census off
  the ORDINATE cosines, never the orbit barycentres (post-#434-R4 those are
  `P_H p`, so the mirror column is exactly zero and `np.sign` reports a
  plausible, flattering, wrong 0 of 4). → `L75c`
- **The SN regression corpus pins NO full-solve `angular_flux`** — the 14 DD
  cases pin `keff` + `scalar_flux`; `2d_octant_equivalence_*` pins angular
  per-SWEEP; the only full-solve angular wall is #448's 32 anchors.
  `material_xs_field()` is a FRESH MINT per call; `geometry_cache_for` fires
  1182× per 1-D eigen solve (intern-absorbed). → `L83k`
- **`-W error::tests.sn.regression._regression_assert.DriftWarning` is a 1-ULP
  wall but NOT absolute on this tree** (19 passed plain, 9 failed / 10 passed
  escalated). Used absolutely it is 9 false reds; used as a DELTA (the drift SET
  and each case's ULP count unchanged) it is exact and free in ~1.6 s. The
  bit-exact cases are all four fixed-source cases (both P1-aniso and the
  windowed 2-D) plus `cyl_2g_3reg_folded_4x6` — the anisotropic and windowed
  paths are the strongest free anchors in the tree. Verify both that the `-W`
  string PARSES and that it bites. → `L58c`, `L77d`, `L83k`
- **MMS in this tree.** An MMS fixed-source is INHERENTLY anisotropic (streaming
  manufactures an ℓ=1 source even for an isotropic trial), so verify a fold's
  MOMENT REACH ≥ the source's anisotropy `L18`. The non-vanishing-at-face family
  LANDED as §4.6 (`build_slab_{,2g_}nonvacuum_mms_case`,
  `build_sphere_nonvacuum_mms_case`, `build_2d_cartesian_ld_stress_mms_case`;
  anisotropic `(A_g + μ_n B_g)/W`) — **do NOT re-derive it; re-route it** `L40a`.
  ⛔⛔ And `vv` Mode 7's "override the simplification bias — high frequency,
  mixed scales" is SCOPED to a SPATIAL-DISCRETIZATION claim: **the strengthening
  axis must be the one the claim lives on.** For a trace claim it is the ANGULAR
  content (`b₀/a₀`); for an angular-CLOSURE fixture it is PARITY, not frequency
  — one EVEN harmonic, then stop, because the τ-independent floor grows faster
  than the signal. → `L40b`, `L48d`
- **⛔⛔ ASK WHAT FIELD MAKES THE SUT'S OWN RESIDUAL ZERO** — a fixture in the
  SUT's kernel cannot rank it, however rich it looks. The shipped curvilinear
  aniso MMS is `A(r) + B(r)η`, affine in the radial cosine, and the M-M closure
  is EXACT on `span{1, μ}` BY DEFINITION of τ, so the flagship angular fixture
  has ZERO closure residual for the scheme it grades. One line of algebra, no
  run. The same check kills the diffusion-limit instrument for ANGULAR claims
  (the diffusion limit's angular content IS `span{1, μ}`) while leaving it sound
  for SPATIAL ones. → `L48a`
- **⛔⛔ When a code path is gated by a parameter's CONGRUENCE CLASS the frozen
  corpus samples one class only** — a whole carve ran on cylinder `n_phi ≡ 2
  (mod 4)` and nowhere else, so `4, 8, 16, 32` READS as a refinement ladder and
  is a single residue; EVERY frozen artifact was blind, including the plan's own
  named canary. Run a counting spy and confirm the changed line EXECUTES before
  crediting any snapshot as an anchor. ⭐ The tree usually already knows: two
  authored comments stated the rule and shipped the activating fixture.
  → `L63a`
- **⭐ Re-run the activation question PER PHASE and PER CLAIM — the answer can
  INVERT inside one campaign**, and the two halves of one step can have DISJOINT
  activating configs (per-cell scheme dispatch: slab 80 / curvilinear 0; closure
  dispatch: slab 0 / curvilinear thousands), so neither geometry family alone is
  an acceptance set. → `L64f`
- **⭐ Production exercises a shared mechanism on a DEGENERATE slice, so the
  general term is never activated — MANUFACTURE the activating case and make the
  load-bearing mutation RED on it and GREEN on production's. That asymmetry IS
  the evidence.** Instances: one seed level makes `pos ≡ 0` `L20`; the S/F arms
  feed ℓ=0 ONLY `L22`; a single-draw probe nulls a two-face law `L32`; a slab is
  the degenerate two-face case for any partner map `L33`; a mint consuming a
  FLAT collection is rank-1 by construction and d=1 hides it `L65d`; a SYNTHETIC
  fixture can null a property the REAL data exercises, making a synthetic-only
  assertion FALSE-RED on production data (pin a cumulative or inequality
  property, never a brittle exact index) `L1`; a SINGLE-REGION mesh starts its
  only region at the origin, so it nulls every defect keyed on the region's
  inner radius and every defect confined to a shell, and a POWER-OF-TWO
  subdivision of a dyadic length makes a float round trip exact (0 of 16 cells
  vs 8 of 23 at 5/7/11) `L89`.
- **⛔ A branch added to DODGE a rank/carrier hazard CREATES the congruence
  blindness** (the new path runs on one carrier kind only), and a second MINT
  SITE hides on the branch where the producer does not exist — a Pattern-2 twin
  one label typo apart. Parametrize the gate over the BRANCH and make the other
  arm assert the OPPOSITE claim. → `L65d`, `L66d`
- **⛔ Before promoting an observed regularity to an assertion, run it on every
  channel the same code path serves** — the (n,2n) and elastic Legendre moments
  decay monotonically, thermal does NOT. → `L76d`
- **⛔ Census which shipped members have an EMPTY channel before letting a new
  length join any `min`** — 2 of 13 isotopes carry no (n,2n), so a two-list
  clamp forces P0 on every water-bearing solve, deleting the ELASTIC P1/P2 (14×
  the effect the campaign existed to add). A control arm must zero the ℓ≥1
  VALUES at the same length, never SHORTEN the list. → `L76a`
- **⛔ A pseudo-inverse round trip is `P_range(G)`, not `id`, and a corpus can
  dodge the null space entirely** — all four SN ledger fixtures carry 0
  tangential (`|Ω·n| = 0`) trace slots while a legal 2-D `product(4,4)` mesh has
  32 of 64, so a naive `== id` gate is blind on the whole corpus AND a false red
  in production. → `L67c`
- **⛔ Before pricing a change to a diagnostic, check whether the tree's pin of
  that diagnostic EXCLUDES the arm the change lives on** — one SI-trajectory pin
  refuses a WINDOWED fixture by construction while the moment metric is read
  ONLY on the windowed arm. → `L80c`
- **⛔ A functional that silently accepts the WRONG SHAPE is more dangerous than
  one that raises, and a 1-group fixture cannot tell them apart** — one consumer
  REFUSED a stored `(ng,)` flux loudly while its sibling ACCEPTED it and
  returned a different number, and at 1g the two spellings agree exactly. ⟹ a
  design handing one datum to two consumers must probe BOTH, and **≥2G is
  required for a SHAPE reason**, independently of the 1-group eigenvalue
  degeneracy. → `L86d`
- **⛔ "Re-point the space" is not plumbing when today's space has NO metric** —
  a `Field.l2` moved 41 %, not ULP. Check `space.inner_product_weights is None`
  before believing any re-point is neutral. → `L62c`
- **⭐ REUSABLE ANCHOR for any ANGULAR-BASIS / moment claim: the infinite medium
  is a Pℓ-ORDER-INVARIANT closed form.** Flat + isotropic ⟹ `φ_ℓ ≡ 0` for ℓ ≥ 1
  ⟹ the anisotropic source is inert ⟹ `k = k_inf` at EVERY truncation order, and
  `derivations.get(...).k_inf` has no solver, quadrature or basis in its chain.
  ⚠ TWO mandatory activation obligations: assert `SigS[1] ≠ 0` IN the test, and
  pose at `scattering_order ≥ 1` — at `L = 0` the folded and parent tables are
  bit-identical. → `L68d`
- **⭐ Subcritical SN slabs are cheap and the anchors already own one** — 2g
  fuel|moderator GL-8 4+4: `L=2.0` refl|vac → `k = 0.435195214`; `L=4.0` →
  `0.907457573` (`1/(1−k) = 10.8`, the strong discriminator); `L=8.0` refl|refl
  → `1.374233987` (the supercritical refusal witness). Three solves 1.17 s;
  three dense 160×160 pencil assemblies 0.23 s. → `L84f`

## 4. Reference, claim layer, and the proactive refutation

- **⭐⭐ A test's CLAIM KIND is the PROVENANCE of its expected value — THEOREM /
  REFERENCE / RECORD — a different axis from `l0`–`l3`** (which grades how GOOD
  the reference is). THEOREM = entailed by a law holding for every admissible
  input (identity, adjointness, involution, conservation, `M − N ≡ A`); red ⟹
  the object violates its own definition and every other claim on that subject
  is VOID. REFERENCE = a structurally-independent external route (`vv`'s three
  pillars); red ⟹ the implementation disagrees with the math *here*. RECORD =
  whatever the code produced on a chosen day; red ⟹ *something changed*, ZERO
  information about which side is right (`numerical-bug-signatures` Signature
  10). ⛔ It cannot be DERIVED: `assert_allclose` appears in 218 files, spelled
  identically for closed form, MMS and frozen baselines. ⭐ The audit it unlocks:
  **every RECORD subject must also carry a THEOREM or REFERENCE test.** Honest
  limit, shipped inside the audit's output: it finds subjects with NO independent
  pin, never a BLIND one — that is mutation's job. → `L55i`
- **RULE: write the (claim-layer, pillar, truth-source) triple per gate BEFORE
  drafting it** (`AGENT.md` §1.5; `vv` §pillars / anti-#5, #6, #7). The standing
  ORPHEUS RESIDUAL: **no mesh-independent transport eigenvalue reference exists
  here** — heterogeneous references are diffusion-based (~0.3 % gap) or
  self-referencing, so the diffusion eigenvalue is a cross-check with an explicit
  tolerance and NEVER a precision target (issue #8). → `L2`
- **⭐⭐ When ORPHEUS has no independent reference, the way past is a different
  CLAIM LAYER, not a weaker gate.** When a solver returns TWO members of one
  object, the reduction identity between them is a FREE L1 gate and the only one
  that can see a defect in the RETURN: `Solution.scalar_flux` is *defined* as
  `∫ Solution.angular_flux dΩ`, needs no external truth, and is a flux-shape
  claim, so the pillar rules hold. It separated by 1.6e6–3.6e6 × its band on 8
  arms at L≥1 while every L=0 control stayed green. → `L78a`
- **Two-anchor template for a pure-refactor carve:** a committed snapshot
  ("didn't move" = bit-id inheritance) is necessary-NOT-sufficient — ULP distance
  cannot tell you the pre-carve value was right — so pair it with a closed-form
  value anchor (`Q/Σ_t`, `k_inf`). Before minting a pre-carve anchor check
  whether the campaign ALREADY froze one (⛔ never RE-capture it), then mint at
  the tier the existing one cannot LOCALIZE: an end-to-end byte capture cannot
  localize the OPERATOR tier, so `A`/`F` against raw `Mixture` arrays is the
  net-new REFERENCE. → `L2`, `L81a`
- **⭐⭐ Two free REFERENCE-class oracles worth reaching for first.** (a) The
  DENSE pencil spectrum: `ρ(A⁻¹F)` from `loss.as_matrix()` and the posed
  `[[F]].as_matrix()` reproduces `solve_sn(...).keff` to 9 significant figures in
  0.23 s — gate it at the solve's own `keff_tol`. (b) Sherman–Morrison for a
  rank-1-multiplying source problem, `(A−F)⁻¹q = A⁻¹q + (A⁻¹χ)(νΣ_f·A⁻¹q)/(1 −
  k_∞)` — `max|Δ| = 0.0`, structurally independent because it never inverts
  `A − F`. → `L84d`, `L84e`
- **RULE (identity-level): the highest-value output of a proactive dispatch is
  REFUTING the plan's optimistic premises with a MEASUREMENT.** Measured false so
  far: "bit-identical"; "clean O(h²) at S16"; "improves on flat at the boundary";
  "this bare-`ndarray` arm is DEAD" (an argument annotated `T` is the strongest
  reason to suspect the `T` arm is LIVE); "N pyright errors clear" (never trust a
  count — assert the residual verbatim); "the same fold applies uniformly across
  N solvers" (the SN/CP/diffusion `keff` DENOMINATORS are different physics); "no
  reported number changes" (a universal over the CALL SITES — 4 of 5 entries were
  bit-exact by object identity and the fifth moved +5.59 %, and the site where
  they disagree is exactly the site the step exists for). → `L10`, `L86a`
- **⭐⭐ A brief's "central risk, ALREADY REALISED in the tree" is a claim to
  audit, and its ENUMERATION is usually short.** When the refutation lands the
  phase collapses from *build a new reference* to *re-route the existing one* —
  also the Pattern-2-correct answer. ⭐ And the ask itself can LAND mid-design:
  run an existence check per promised DELIVERABLE, not only per named symbol.
  → `L40a`, `L43a`, `L43j`, `L39`
- **Measure the proposed ACCEPTANCE CRITERION as a probe before any gate is
  written** — an AC shaped "changing X must not touch Y" is usually already true
  BY SIGNATURE, so it is unfalsifiable from the first commit and a falsifier
  check PASSES on it (`vv` Mode 8(3)). Gate the SIGNATURE; demote the value row
  to a regression floor. → `L24`
- **⛔⛔ A brief's headline NUMBER carries an unstated REGIME — reproduce it
  before designing to it, and say so if it only reproduces off the production
  path** (a reported `min ψ̂ ≈ −77` reproduced only with a RANDOM ψ and a ZERO
  seed; on the production value path the same fixture gives **+0.13**). Then pin
  the MECHANISM, not the observation: the mechanism was solve-free, a pure
  function of the chart, with a closed-form independent reference explaining BOTH
  regimes. → `L47b`, `L47c`
- **⛔⛔ A brief saying a relocated computation uses "plain / flat / simple"
  arithmetic has named TWO conventions — enumerate the candidate spellings and
  MEASURE the spread before writing the pin; the spread IS the pin's
  discriminating power and its tolerance.** The dangerous rival is the CONVENIENT
  one (it deletes a helper). ⚠ A "harmless" rival can be phase-scoped, so make
  the later phase's legitimate red a blocking ruling, not a surprise. → `L58a`
- **⭐ A "the fall-through lands on X after we delete Y" claim is a RUNNABLE
  experiment, pre-carve** — invoke the base implementation the carve will expose
  on the inputs the guard must still refuse. ⚠ Assert the DISCRIMINATOR in-test:
  `a.space == b.space` can be True (`FunctionSpace.__eq__` is `(name, shape)`)
  while `a.space is b.space` is False, so a row omitting that precondition
  degrades silently. → `L58d`
- **⛔⛔ "Solve for X on interval I" is a CONSTRUCTION only if the root is UNIQUE
  there — run a sign-change scan and report the root COUNT before accepting the
  design** (a briefed root-find had TWO roots at 4 of 9 orders, so `brentq`
  raises before it starts). Two roots ⟹ the plan owes a SELECTION RULE with a
  two-legged gate: (a) the shipped value IS the selected root, (b) the discarded
  root is exhibited and measured bad — (b) is what makes the rule a reason
  instead of a coincidence, and the only attributable leg. → `L42a`
- **⭐ When a plan's blocker is "is this small number real or is it arithmetic?",
  answer with arbitrary precision instead of escalating** — one mpmath probe at
  50/60 dps closed a blocking user ruling. Keep the reasoning ✅ ANSWERED, not
  deleted (`plan-authoring` §3), and pin the MARGIN VALUE beside the sign so a
  conditioning regression reds first. → `L42e`
- **⛔⛔ AN ADJUDICATING INSTRUMENT (one that RANKS designs) is a different object
  from a gate**; its four checks are now `vv` anti-#24, promoted from here. What
  anti-#24 does NOT carry: **a CONTINUOUS homotopy beats "rank agreement"** —
  require MONOTONICITY in `w` along `blend(w) = (1−w)A + wB`, five solves,
  falsified by one non-monotone triple; STRATIFY the ensemble (NEAR = a 2 %
  jitter · MID = the real rival · FAR = garbage) and require the threshold on
  **NEAR∪MID alone**, since a ρ over all three is dominated by the garbage split;
  the ensemble MUST contain the pair inside the stabiliser you fear. Declare each
  instrument CONSTRAINT / RANKER / DIAGNOSTIC in its own docstring — the
  graveyard died of silent promotion, and 4 of 6 dead instruments died at the
  `<1 s` solve-free pre-flight. → `L48a`, `L48b`, `L48f`
- **⭐⭐ "This comparison is BELOW MY RESOLUTION, and here is the number" is a
  first-class deliverable, not a failure.** The best fixture and functional
  resolved garbage 17–40× and a 2 % jitter 2–4× — and NOT the two candidates the
  campaign actually argued about, which makes "decide on constraints + the
  primary source" the sound route rather than a fallback. ⛔ Related:
  **closure-EXACT is not accuracy-optimal**, so every closure-residual instrument
  is a DIAGNOSTIC. → `L48c`, `L48g`
- **⭐ The keystone's ORACLE choice decides whether it catches anything — same
  assertion shape, 8 orders of sensitivity apart.** For any "the answer satisfies
  the declared condition" gate ask **which side is the thing under test**; if the
  answer is "both", it is not a gate. A rewire's demotion test is the same
  question: **is the retired symbol a SOURCE of the expected value or a FORWARDER
  of it?** → `L40e`, `L39`
- **⛔ A re-pose can INVERT a migration gate's SENSITIVITY partition — the
  inherited `[M]` characterisation dies by being FIXED, not refuted.** The old
  anti-claim arm becomes a must-RED arm and a brand-new must-stay-GREEN arm
  appears (the *un-wiring* proof) that could not be stated before. Run both at
  BOTH HEADs and put the 2×2 in the docstring. → `L61c`
- **⛔ A corpus paragraph can carry an honest `[M]` whose LOAD-BEARING half is
  false, because its experiment varies two things at once** (anti-#17(a)'s
  granularity trap at the doc tier) — and a carve can make the claim TRUE and its
  mechanism clause present-tense-FALSE at the same time. → `L61d`
- **⭐ Before excluding a field from an identity key on DOCTRINAL grounds, check
  whether the exclusion is also MANDATORY** — the stronger, more durable gate.
  Including a `Quadrature`/`DiscreteMeasure` makes `__eq__` RAISE and `hash`
  RAISE, not merely disagree; pin the REASON with `pytest.raises` legs on the
  generator TYPES. → `L65b`
- **⚠ Name which half of a comparison is REAL.** Two "mint vs literal" gates read
  the SAME array object on both sides, so they pin THREADING (label, shape
  spelling, `kind`, wiring) and never the values; the honest digest gate rebuilds
  the pre-change literal space IN THE TEST. And a bit-identity row comparing two
  BINDINGS OF THE SAME CLASS cannot see a defect inside that class (anti-#22's
  third manifestation: no shared object, no caller relation, just a shared
  implementation) — for any "two mints agree" row, mutate what makes them DIFFER.
  → `L65g`, `L83g`
- **⭐⭐ Measure an ADJOINT/metric objection on the RANGE OF THE PRODUCER, not on
  `randn(space.shape)` — a claim about inputs the producer cannot emit is not a
  claim.** A recorded "moves 10 of 33 rows" was 5 of 33 on random draws and 0 of
  33 on a covariant moment `φ = Mψ`, because on a folded rule the σ-odd harmonics
  are identically zero at every node so `G⁺` projects them out. For any `.H`
  claim on a producer's CODOMAIN, the fixture is `producer(x)`. → `L80a`
- **⛔ An inherited `[M]` PERCENTAGE with no statistic is unreproducible —
  replace it with the DRAW-FREE one rather than hunting for the original.** For a
  diagonal-metric swap the honest statistic is the per-element ratio `|p_i/g_i −
  1|`; an L2 residual is draw-dependent. → `L80b`
- **⭐ A registry field asserting a PHYSICS claim can be gated at the SOLVER tier
  for ~1 s, and that is the most a gate can say** — solve a deliberately
  asymmetric fixed source and compare ψ at ordinate n with ψ at the ordinate g
  maps it to, one positive leg per element IN the group and one NEGATIVE leg per
  element outside it. ⭐ And do NOT gate a table relation you cannot DERIVE:
  record it as an observation WITH its denominator, nowhere as an assertion.
  → `L75c`, `L75`
- **⭐ When a census says a whole FAMILY is blind, look for an existing REGISTRY
  case before designing machinery** — one shipped case closed a 0-of-113 gap in
  0.005 s. And a stochastic method's "too slow to gate" is usually a statement
  about the PRECISION target, not the catcher: a 0.9 s MC run read 0.47 σ honest
  and 17 σ mutated where the only `slow`-marked catcher is deselected by the
  canonical `-m "not slow"` (`vv` anti-#36). → `L76g`, `L76h`
- **⭐ A declared ONE-SIDEDNESS needs its own battery arm — GREEN on the blind
  gate and RED on its two-sided partner.** The physics bound `|Σ_ℓ| ≤ Σ_0`
  catches an inflation and is blind to a deflation; the two-sided catcher is a
  RATIO-INVARIANCE row, and the arm proving the pair is `scale**ℓ`. ⚠ Threshold
  `1 + 1e-9`, never tighter. → `L76c`

## 5. Tolerance is a claim — choose it per law, from measurement

- **RULE: bit-exactness is EARNED PER LAW; measure before choosing the
  assertion.** On ONE type: identity 500/500 bit-exact; associativity 500/500 on
  signed permutations, 0/500 on general rotations; `g∘g⁻¹` 0/500. A uniform
  choice is a false red or a thrown-away gate. → `L35h`
- **A law's BIT-EXACTNESS can be ARM-DEPENDENT and FACE-DEPENDENT — measure per
  arm and per face before writing `array_equal`.** One splitting residual is
  exactly 0 on every SEEDLESS arm and 3.6e-15…2.8e-14 on the CARRYING arm (the
  grid re-associates); a forward READING can be `array_equal` 200/200 while the
  forward APPLYING on the composite is 0/200 at ≤1 nulp. Name the draw-stable
  statistic for the non-exact arm. → `L77c`, `L83a`, `L84k`
- **State the law in the direction that IS a float theorem, and normalise a
  residual that scales with its input.** `on_points − on_directions == t` is NOT
  exact (`fl(a+t) − a ≠ t`); `on_points == on_directions + t` IS bit-exact
  6000/6000 because it recomputes the same expression — and is the stronger
  assertion. An ABSOLUTE `atol` for a residual scaling `O(ops × ‖t‖ × eps)` reds
  on correct code for large draws and is too loose for small: normalise by
  `max(1, ‖desired‖_∞)`. → `L35j`
- **⛔ A DIFFERENCE of two members of an affine operator family is NOT
  bit-identical — the MATRIX form is.** `at(σ+τ) − at(σ) == −τ·M` differs on 144
  of 200 draws (catastrophic cancellation, not an identity) while
  `at(σ).as_matrix() == A − σ·F` is `array_equal`. Same family: `Σ(Aψ) − Σ(Sψ)`
  is NOT `Σ((A−S)ψ)`, so a re-signature folding two operands into one changes an
  estimator's arithmetic — and it can hold EXACTLY on a synthetic fixture whose
  hand values are representable, so the committed row survives the carve for a
  reason that does not generalise. → `L84b`, `L85f`
- **Re-derive every tolerance from structure; retire inherited `nulp` folklore.**
  Gathers and α-folds are reduction-depth 0 ⟹ `array_equal` (a tolerance there
  would admit the bug); an `n`-term positive-summand contraction vs a `tensordot`
  is `κ=1` ⟹ `|Γ₊|·ε`, and the probe being non-negative is WHY `κ=1` — say so.
  For a "hand it the constant" move the realistic defect is the CLEANER algebraic
  spelling, 1–2 ULP: `array_equal`, since any tolerance ≥ 1e-15 is a non-catcher.
  → `L32`, `L63g`
- **Regression-snapshot tolerance is the CLAIM, not a magic floor**: iterative →
  `SAFETY(10) × conv_tol` read OFF the run config (the SoT shared by generator
  and test); direct → `nulp(reduction_depth)`; bit-identity by `-W
  error::DriftWarning` LAYERED on top. Corollary: an ITERATED end-to-end snapshot
  CANNOT be the bit-identity gate for a zero-numerical-change refactor —
  committed iterated snapshots already drift 1000s–100000s ULP from cross-run FP
  jitter; descend to a single-step DIRECT snapshot on a fixed-seed random
  heterogeneous ≥2G ψ with non-zero inflow. Recipe:
  `feedback_regression_tolerance_design.md`. → `L7`
- **⛔⛔ For a ROOT-FIND gate the tolerance is `noise / slope`, and the slope can
  collapse 4 orders across ONE parameter family** — so a single rtol is a false
  red at one end and a dead gate at the other (1.0 ULP at S4 → 40 653 ULP at
  S18). Derive it (`Δx ≈ evaluation noise / |f'|`), tabulate PER ROW ×10
  decade-rounded, put the arbitrary-precision value in the literal, and STATE
  what the floor leaves ungated. → `L42b`
- **⛔ A flat `atol` is wrong in BOTH directions at once — derive it from what the
  quantity DIVIDES BY, and note that two quantities in one seam can need two
  different laws.** One τ row carrying `atol=1e-13` was ~450× too loose at N=8
  AND a false red at the order its own docstring predicted. Derived: sphere τ
  divides an `O(ε)` edge discrepancy by the cell width ⟹ `16·ε/w_min`; cylinder τ
  inherits `cot`'s conditioning ⟹ `40·M·ε`; while the PARTITION the same τ reads
  agrees at a flat ≤1.5 ULP. Same for a negative control's FLOOR: one convention
  gap SHRINKS like `M⁻²` in edge space and GROWS in τ space. → `L47f`, `L47h`
- **⛔ When a guard compares two independently-accumulated floats the tolerance is
  a MEASUREMENT over the constructible population PER SUB-FAMILY**, never a
  judgement about whether the construction "should" be exact — 0 ULP on slab and
  every `uniform` mesh, 1 ULP on CYL/SPH `equal-volume` (a `sqrt`/`cbrt`
  round-trip), so `==` is a *latent* false red. Ship a derived band WITH its
  discrimination margin, and put the arm that PROVES it in the battery. → `L60b`
- **⛔ "Bit-identical at the degenerate fixture" is usually 1 ULP — asserting
  `array_equal` on it reds your OWN control** (`np.cos(np.pi/2) = 6.12e-17`, not
  0). Assert "15 orders below the signal", never "the bits match"; and such a
  blindness CONTROL still earns teeth — it must red when one convention moves and
  stay green when production BECOMES the other one. → `L47g`
- **Gating a MORE-accurate implementation against the less-accurate one it
  replaces gates it against the error it exists to remove.** State the criterion
  against the arbitrary-precision value, and honestly ("within 0.57 ulp
  everywhere", never "always closer"). → `L34c`
- **⭐ An instrument's ORDER is a measured choice, and MORE CAN BE WORSE** — an
  orbit-circle trapezoid mean is exact for `n ≥ 3`, so its residual GROWS from
  2.2e-16 at n=8 to 2.6e-14 at n=1024. Ship the small `n` and SAY in the
  docstring that raising it degrades the gate, or a later session "strengthens"
  it into a false red. → `L73h`
- **⭐ A random-draw separation statistic usually has an EXACT draw-free
  replacement — ask whether it is a Rayleigh quotient.** A committed floor pinned
  a SEED; the same statistic ranges 0.23…2.00 over 400 draws on the very frame it
  gates, and its exact range is one `eigvalsh`. → `L69d`
- **⚠ An asymmetric-morphism law needs its activation PRECONDITION ASSERTED** — a
  condensation pair discriminates against its three wrong pairings only while
  every coarse group holds ≥ 2 fine groups; at one fine per coarse two of three
  controls go silent. → `L67h`

## 6. Carve archetypes — where the load-bearing gate lives, by carve shape

Reference material, not a per-dispatch rule: moved to `lessons_archive.md`
§L87 on 2026-09-21 (60+ shapes indexed by their bold name — relocation,
un-weld, type-collapse, admission/refusal, rename-a-field-and-move-its-
meaning, kernel/datum-mint, ...). Skim §L87's bold names when a dispatch's
carve matches one; open only that row.

**Meta-rule (kept hot — it is the one habit, not the lookup): the keystone is
decided by whether the carve INHERITS a verified predecessor.** Wrapping /
re-expressing something verified ⟹ the keystone is bit-id INHERITANCE
(necessary-NOT-sufficient, always paired with an independent value anchor).
Nothing to inherit ⟹ the keystone must be structurally independent. Before
accepting any BIT-IDENTITY acceptance line ask which REDUCTIONS the change
reorders: zero ⟹ `array_equal` is honest; any ⟹ the line is arithmetically
IMPOSSIBLE and must be re-scoped to a permutation reordering no addition.
→ `L37`, `L87`

## 7. Snapshots, generators, and exactness

- **RULE: a snapshot generator that calls production and freezes its output is
  SELF-REFERENTIAL** — it says `production == a recording of production`, detects
  change and certifies nothing. INVERT it: compute the reference from the law's
  EQUATION (never by transcribing the implementation) and freeze THAT;
  precondition, the expression must be TOTAL. Then the FROZEN FILE is the only
  thing between a wrong expression and a green gate — make that structural: an AST
  gate asserting the generator imports nothing from the realization layer, the
  harness pulls only the case registry, artefacts on disk == registered cases.
  → `L32`
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
  ambiguity converge away? (4) SPLIT the commits, the ordering ruling first,
  alone. And **"the level is sorted by η" is the wrong functional** — sortedness
  is invariant under permuting equal elements. Gate the full INDEX TUPLE against
  an independently-constructed one, plus a `kind=`-invariance row (quicksort /
  stable / heapsort / mergesort must agree bit-identically), the operational proof
  the key is injective. → `L34`
- **A nearest-neighbour partner search is not automatically a bug — measure the
  MARGIN first** (separation 5.0e-3 vs a 1e-16 perturbation REFUTED the claim);
  keep the search rather than minting an index-formula twin path. And `ref[ref] ==
  id` passes on a residual-0.94 garbage map: any self-inverse pairing satisfies an
  involution law, and the only functional outside that stabiliser is the RESIDUAL.
  → `L34b`, `L33`
- **⛔⛔ A principled re-read can sit INSIDE the hard tolerance band and OUTSIDE
  bit identity — name the DRIFT TRIPWIRE's new population, not just the gate's
  headroom.** Deriving the eigen `scalar_flux` as `∫ψ dΩ` was `array_equal` 0 of
  16 cases against a band with 135× headroom (every pin holds), while the
  escalated `-W error::DriftWarning` run goes 2 → ~28 reds and STAYS there. A
  bit-identity claim in this tree faces TWO gates; say which one you mean. → `L86f`
- **⭐ Probe the algebra before choosing a tolerance — gate-ready bit-exact laws
  are FOUND by probing, not assumed** (`R∘E = id`, `E∘R` idempotent, `R.H == Σw·E`,
  `analyse(ℓ=0) ≡ integrate_angular` — all `np.array_equal`, four needing no
  tolerance). And when a new type must reproduce an existing table bit-exactly,
  diff the existing producer's BRANCHES, not its name: no single library routine
  reproduces a hand-branched table, and a spot-check at one ℓ certifies the wrong
  spelling. → `L62d`, `L69a`
- **⭐ When a carve only RELOCATES where a measure is stored, simulate it with
  `dataclasses.replace` pre-carve and claim `array_equal` — but MEASURE it.**
  `dataclasses.replace` round-trips a frozen canonicalizing dataclass
  bit-identically (bytes, read-only flag, `eq`/`hash`, kw-only subclass fields,
  idempotent canonicalization), so a `replace`-based field upgrade is
  `coding-elegance` Pattern 4∩2 safe and needs no hand-written constructor.
  → `L65h`, `L80d`
- **⛔ A derived space NAME is a future landmine when axis identity is
  per-SUBCLASS** — pin axis CONTENT and relative identity, never the name literal.
  ⭐ And an identity derived from a RENDERING is only as injective as the
  renderer: ndarray `repr` TRUNCATES with `...`, so two distinct long weight
  vectors give the SAME name. Derive from `.tobytes()` through a digest, state the
  float caveats (`-0.0` vs `+0.0`; `nan`), and give the injectivity gate at least
  one pair whose SHAPES are identical. → `L59e`, `L62f`
- **⛔ "Prove it does not allocate" cannot be gated by asking a densifier to
  `MemoryError`** — a 550 GB outer product is OOM-KILLED (exit 137), which fails
  the RUN, not the test. Size the gate for SEPARATION (asserted on reachable
  `ndarray.nbytes`), add an EXACT structural leg, and add a BEHAVIOURAL leg —
  "never densify" implemented by DROPPING the metric passes the first two.
  → `L59c`

## 8. Verifying a pure-math PRIMITIVE (a group / algebra type)

`vv-testing` § "A type that embodies a mathematical concept ships the test of its
defining laws" is the obligation. Reference material for the fourteen further
gate shapes (homomorphism laws, involution/affine blindness, the G-fixed
centroid, bijectivity vs match window, the coset-search inverse-direction
theorem, the stabiliser-maximality gate, ...) moved to `lessons_archive.md`
§L88 on 2026-09-21 — open it when the dispatch is a group/algebra-type gate.

**The one fact worth keeping hot: the pillars differ from a solver's — no MMS
row, no semi-analytical row.** Every row is closed-form; the structurally
independent grounds are SymPy under an EXPLICIT unit parameterisation
(imposing `Σnᵢ²=1` by `subs` after expansion does NOT fire), an external
implementation with a DIFFERENT ALGORITHM (quaternion vs Rodrigues), the Lie
definition `expm(θ(vuᵀ−uvᵀ))` (dimension-generic, ~4e-14 ⟹ gate at 1e-12),
published tables, and EXACT INTEGER arithmetic — the last needs no reference
at all and is the strongest class available. → `L35a`, `L88`

## 9. Pointers

- **Characterization vs guarantee:** GUARANTEE tests carry `verifies(...)` and
  assert what IS correct; CHARACTERIZATION tests carry NO `verifies(...)` and
  bound a limitation ONE-SIDED (no upper bound, so a future fix keeps them green).
  To pin a floor a fix claims to remove, measure the floor's SCALING with the
  OTHER axis — `err(S32) < err(S16)/2` is falsifiable where "the floor is gone" is
  not. An out-of-scope defect gets a POSITIVE assert-the-defect gate with a loud
  message, NOT a non-strict or imperative xfail; a `strict=True` MARKER is the
  preferred spelling, since its XPASS is a FAILURE and it retires itself. → `L5`,
  `L16`, `L45`
- **Mode-10 sub-floor terms:** producer-threading at machine precision + a
  consumed-flip ≫ tol + a no-op control; where NO isolating regime exists the
  ABSENCE of a value-improvement leg is the CORRECT signature (`vv` Mode 10).
  → `L6`
- **Two instrument vocabularies, do not mix them.** THEOREM / REFERENCE / RECORD
  (§4) is a TEST's claim kind. CONSTRAINT / RANKER / DIAGNOSTIC —
  `docs/development/evidence/lessons.md` § "L51 instrument type" — is for a
  DESIGN-RANKING instrument, and is the wrong vocabulary for a test (every
  collected pytest test is a CONSTRAINT by construction, so that partition has one
  non-empty cell).
- **V&V tagging idioms** (`foundation` must not carry `verifies`; module
  `pytestmark`; `slow` on params not functions) → `feedback_vv_tagging.md` (`L9`).
  **Cross-method agreement infra** (reuse the registry schema; `max(tol_a,
  tol_b)`; tag L1 not L4; truth values MUST trace to primary citations) →
  `feedback_cross_method_protocol.md` (`L8`). **Per-carve RECIPES** → `MEMORY.md`
  §3. **Measured scope costs and the serial-host gate** → `MEMORY.md` §2 and the
  main memory's `reference_test_execution_env`.
