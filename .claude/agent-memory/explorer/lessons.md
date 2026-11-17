# Explorer — Lessons (the hot digest)

Behavioural corrections only: *what mistake did I make exploring, and what
changed how the NEXT exploration is run?* Each entry is an imperative, its
tell, and a pointer. The war story behind an entry numbered L-001 to L-045 is
cold, in `lessons_archive.md` under the same L-number: open a section only
when a pointer needs checking. L-numbers are stable and never reused. An ID
absent from §2 is either wholly carried (by the definition, by a preloaded
skill, or by a generated `brief:` line) or survives as one clause of §3.

What reaches this agent: its definition, the three preloaded skills
(`nexus-exploring`, `nexus-guide`, `retirement-audit`), and the generated
`brief:` lines of seven rules when the brief's "Rules that apply to you" list
carries them whole (`articulation`, `code-search`, `coding-standards`,
`instrument-doctrine`, `plan-authoring`, `process-discipline`, `workflows`).
A rule's BODY does not reach it, so a clause that lives only in a rule body is
not a carrier here; such an entry is marked `[rule body]` and stays.

Two IDs are irregular. Archive §L-046 is the graph-degree lesson, which the
definition §4 carries ("A hub is not a template"); the entry L-046 below is a
different lesson written later under the same number, with no archive section.
L-047 to L-050 have no archive section either: the entry is the whole lesson.

## 1. The meta-lesson that no carrier holds

- **M-3 A census is an instrument: split its populations before counting,
  and control its zero.** The `instrument-doctrine` brief line carries X1 and
  X2 (predicate, tree, exclusions, `k of N`, a positive control per shape, an
  AST pass for membership). The splits this corpus keeps teaching: executable
  call vs prose citation vs duck-typed surrogate (`def <method>` outside the
  owning module) vs the `setattr`/`patch.object` string; receiver-resolved vs
  bare token (a method on one class and a field on another share a spelling,
  so collect the receiver per hit); xref vs promissory claim vs convention
  statement; anchored vs unanchored (the English-word substring, the
  same-package homonym); the AST node-type viewport (`Assign` misses
  `AnnAssign`, `FunctionDef` misses `AsyncFunctionDef`, `Name` misses
  `Attribute`: a zero over a known-populated file indicts the predicate
  first); two counters at different depths so "built but never reached" is a
  visible row. A census that disagrees with a brief by ONE is almost always a
  prose line or a surrogate. Instances: L-009, L-030, L-043, L-049, L-050
  here; archive §L-017, §L-027, §L-037, §L-038, §L-039, §L-040, §L-045.

M-1, M-2, M-4, M-5 and M-6 retired into the definition on 2026-09-22 (§3, §3,
§2 last paragraph, §2 last bullet, and Return, in that order). A pointer
below names the definition section, never a retired M-number.

## 2. The lessons

- **L-003** A transient line map goes into a topic file stamped with its
  HEAD, and is marked re-derivable; the durable shape leads (definition,
  Return and Writes). Archive §L-003.
- **L-006** "Collapse these N shape probes into one predicate": classify each
  as boolean PRESENCE (swap to the typed predicate) or integer WIDTH (an
  honest count feeding allocation: keep), and check whether the probe site
  even HAS the typed object in scope; report "(B) small plumbing", not "(A)
  clean swap", when it does not.
- **L-008** `[rule body]` zsh: an unquoted word starting with `=`
  (`echo ===`) is a command lookup that aborts the WHOLE compound command,
  silently losing every grep sequenced after it. Quote separators
  (`printf 'NAME\n'`, `echo "---"`). The clause is in `code-search`'s body
  and not in its brief line.
- **L-009** A dataclass-FIELD rename is a grep problem (definition §2: no
  field node). Grep `\.<field>\b|<field>=|<field>:` plus `replace(`,
  `getattr`, `asdict`. Before any replace strategy, check whether the old
  token is a substring of an English word (`gains` ⊂ `against`, `loss` ⊂
  `lossless`, `role` ⊂ `payroll`): grep with and without `\b` and hand the
  delta to the implementer as a hazard line. → M-3.
- **L-013** "What breaks if this numeric primitive changes?": swap it (a
  `pytest_configure` plugin) and run the consuming suites; about 200 grep
  candidates became 2 measured items. Classify by "is the RHS FROZEN"
  (`.npz`, a literal, a hash), never "is the comparison exact": same-process
  route-equivalence is immune. Patch every re-export AND every captured
  function object (`object.__setattr__` on the registry field), and print a
  confirmation line. Run the sibling suites too (the only moving snapshot
  was reached through `_generate_snapshots.CASES`). A guard's FIXTURE
  ENUMERATION is where vacuity hides (a gate enumerated GL and Lebedev, not
  `product`); check whether the "new" hazard already exists on the
  already-exact sibling.
- **L-014** Adjudicating an algorithm against the literature: read the
  paragraph that PRODUCES the equation (Hébert defines `α ≡ 𝒲·η` two lines
  above the recursion, which turns "convention?" into a decidable closed-form
  check) and grep the sidecar for the INDEX-DOMAIN prose (`interval`,
  `normalized to`, `for m = 1 … M`). A cumulative recursion telescopes under
  every permutation, so its closure, sum and step gates are ordering-blind.
  Find the closed form before endorsing an MMS run.
- **L-015** "This DOF is redundant, fold it": enumerate the FUNCTIONALS of the
  state (moments, currents, leakage, inner products, the adjoint metric) and
  classify each integrand's parity under the group; odd-parity ones are OUT
  OF THE SPACE, not inaccurate (a fold that was 5e-16 on every even moment
  turned the odd one into +2.94). Split fold-the-ALGORITHM (lift partners
  back; symmetry by construction; no memory win) from fold-the-STATE before
  scoping. When both candidates satisfy the issue's criterion, find the
  structural predicate the tree already keys on.
- **L-018** A hand-written lookup table behind a computed fast path has dead
  rows: per row, "which branch answers this?" (2 of 5 were live). A bare
  `return False` fallthrough gives a new tag a wrong-but-silent answer. A tag
  with two dispatch branches is invisible on a fixture both accept: find or
  build the input that BREAKS the property (`product(4,3)` is closed under
  σ_z and not σ_x).
- **L-019** Hunting a hidden transformation: grep the chart-defining
  ASSIGNMENTS (`cos_theta =`, `= arctan2(`, `nodes[:, k]`), reconstruct the
  implied frame matrix, and ask the checker whether it can NAME it
  (`cos_theta = mu_x` is a 120° `O_h` element no tag expresses). COUNT a
  partition's parts before believing its name (`octants` gave 26 on
  Lebedev(17)). The tolerance-family census is free (three idle epsilons,
  one justified by a symbol that existed nowhere). Evaluate a "for the
  degenerate case X vanishes" docstring ON the degenerate case.
- **L-022** "What REMAINS of issue X": grep the campaign's own interim
  vocabulary ("today still", "remaining half", "until then", "not yet", the
  `.. caution::` directives in touched modules) and re-verify each hit; the
  honest mid-flight notes are what nobody returns to (four falsified in one
  pass). "Did the gate land?" reads the gate's FIXTURE ENUMERATION, not its
  assertion.
- **L-029** A "build the capability vs refuse the input" fork: price the
  capability's OUTPUT against the incumbent at equal cost (a marginalised LS
  rule lost seven orders to `gauss_legendre`). "Does the primitive exist?"
  is per TIER: read the candidate's MATCHING and precondition code, not its
  docstring's verb. A guard's QUANTITY may be annihilated by the input
  class's symmetry (the first moment is 0 on every `O_h` rule). Sweep every
  shipped rule through the real constructor and read the traceback FRAME.
- **L-030** Before instrumenting a seam, `grep -n "return <Type>(\|return
  _<tail>("` and patch EVERY site: a "the ONE construction site" docstring
  was false (10 became 24 affected tests, silently). Build the census with
  two counters at different depths so "config built, seam never reached" is
  a row you must explain. → M-3.
- **L-032** Equation-implementer hunt: the `vv-status rationale` naming the
  answer usually sits on a NEIGHBOURING label (the test's `pytestmark`
  comment, a sibling's sentinel, another page): grep the LABEL across `docs/`
  and `tests/`, not ±20 lines. A NEGATIVE sentinel ("carries no vv-status
  because …") is evidence a `vv-status` grep misses. Before
  `NOTHING:identity`, ask what object the LHS is and grep for a routine
  returning it (`balance_residual`). `equation_labels=(...)` on case
  dataclasses are existing declarations. A dataclass field resolves
  `py:attribute:`, illegal as an `implements` source: escalate to the class.
  KIND prior: a page of ALGORITHMIC RULES is nearly all declarable, an
  ALGEBRAIC-LAW page about half (`{identity, law, canonical-form}` → NONE;
  `{typing-rule, definition}` → look for a declaration site); the measured
  shares are in the `nexus82_*` topic files and archive §L-032.
- **L-033** The CODE may already declare the label (`grep -rn ":label:
  <name>" orpheus/`; `grep -rn "(Eq. <name>)" orpheus/derivations/`), and a
  page may say so. The `@pytest.mark.verifies("<label>")` claimant's BODY
  names the SUT. The principled NOTHING's tell is *independent* / "NOT the
  production" (declaring it would make the gate a self-comparison). An
  equation asserting `A ≡ B` can be HALF-retired: existence-check both
  sides. A `w`-generic primitive legitimately implements several labels:
  say so and let the declarer rule on breadth. Tier the answer (arithmetic /
  generic-at-a-constant / factory / symbolic algebra-of-record).
- **L-034** "Doctrine A replaces B" needs THREE inventories: B's
  implementation; A's existing fragments (the flag nobody reads, the
  refusal, the sibling family's battery, the normalisation; these flip
  BUILD to UNIFY or CONSUME); B's own concessions to A.
- **L-042** Dataclass introspection lies: a `frozen=True, eq=True` record with
  an array field has a generated `__hash__` that RAISES and an `__eq__` that
  raises on two or more elements. Only `hash(a)`, `a == b`, `a == a` at
  runtime answer, on an `ng ≥ 2` fixture (`ng=1` is a false green). Report
  the introspected flags only as "what introspection says".
- **L-043** A kwarg's move census is an AST keyword census split BY RECEIVER
  (148 sites against 183 grep lines; 54 `La13511Case(scattering_order=…)`
  were a record-FIELD homonym). Trace the VALUE per entry to where it is last
  transformed: a knob "lives" where it is last transformed, and per-entry
  disagreement makes the move a twin-path repair. Grep the string class
  (`"X"` in dicts, `cfg.get`) separately. → M-3.
- **L-044** A "bit-identical" verdict: state each wall's PREDICATE class
  (band vs `nulp`/`array_equal`) beside it. ULP-probe the two spellings at
  the operator level on a random operand (they were 2 to 5 ULP apart) and
  name the pre-carve `nulp` anchor that must be captured. Spy pitfalls: patch
  the CALLER module's binding (the definer's binding read 0), key the spy on
  BUILDS not lookups (207 against 1), assert the exact count.
- **L-046** A coverage map built while a peer writes the successor gates
  goes stale DURING the session (line numbers shifted about 40 lines; a new
  untracked gate file appeared after the census). Take `git status --short
  tests/` at open AND at close, run `git ls-files --error-unmatch` on every
  file you call landed, re-derive every cited `file:line` from its anchor
  last (a script over (file, anchor) pairs), and mark each successor HEAD vs
  working tree. The definition's "re-run the empty searches" sees content,
  not a NEW file. → M-3.
- **L-047** "Would hoisting this lazy import create a runtime cycle?": answer
  by an in-memory hoist probe (a meta-path loader that inserts the import
  statement at the top of ONE file's source; one fresh interpreter per import
  order), never from the module graph alone. The group-level graph
  underestimates: a package `__init__` executes its siblings, so 5 of 7
  "needs a lazy hop" typing edges still failed at load time. Tell: the loader
  must override `get_code`, or the cached `.pyc` skips `source_to_code` and
  every probe passes silently; print a rewritten/NOT-REWRITTEN flag as the
  control. Script (untracked): `scratch/numerics_untangling/hoist.py`.
- **L-048** `[rule body]` `git grep -E` here honours neither `\b` (the
  `code-search` body says so) nor `\w`: `class \w*Space` reads 0 where
  `class [A-Za-z_]*Space` reads the classes (`[M]` 2026-11-16, re-run).
  Use `git grep -P`, which honours both, or ugrep with an explicit class,
  and run a sibling pattern that MUST hit as the control in the same call.
- **L-049** A shape-membership census over a space's verbs must separate
  production traffic from GATES that pin a convention: an attributed run
  (record `PYTEST_CURRENT_TEST`) showed every non-member call came from 7
  legitimate tests and 1 mutant, which is the evidence a "refuse at the
  verb" design needs; the unattributed count alone read as "non-members
  exist". Plugin (untracked):
  `scratch/layer_dag/ipf_fix/explorer_member_plugin2.py`. → M-3.
- **L-050** A runtime per-class verdict from SAMPLES is only as good as the
  sample key: a cap of 12 samples per (class, verb) was spent on the angular
  frames and hid the one instance of `FrameAnalysis` that differed. Key
  samples by INSTANCE (parts, shapes, space), or re-run the suspect test
  files with fresh samples, before calling a class uniform. Exclude rows
  whose call RAISED from any verdict: refusal gates feed misfits on purpose.
  A signature row's `first` test names one test, not the call population;
  to attribute production vs fixture, record the 3-frame caller chain
  (`sys._getframe`). Plugin (untracked):
  `scratch/true_ends/census/plug/ends_spy.py`. → M-3.

## 3. Residual checks of retired lessons

Each clause is a `check:` that its retired lesson held and that no carrier
states (restored 2026-11-16 after a re-count of the 2026-09-22 retirement);
the war story is the archive section of the same number.

Sharpening the definition §2 (a behavioural question is answered by a run):

- **L-010** "X is the complement of Y, hence X == Z": find the PREDICATE that
  defines the split and check that it is exhaustive. A strict `<`/`>` pair
  with an epsilon never is, and whether the third bucket is populated depends
  on the data.
- **L-041** An identity predicate spelled with `is` on an `Optional`
  constituent is vacuously true where the constituent is absent: run it on
  the `None` arm AND on a rebuild from the same generating data.
- **L-024** A solver's nesting shape is a fact per PUBLIC ENTRY, never per
  module: enumerate the entries and trace each by running (a `sys.setprofile`
  probe filtered to the files watched), then read each count against its
  budget. A status derived from a LENGTH is checked for injectivity.
- **L-025** "What is in scope at call site P, could P compute Q?": a
  frame-locals spy with the ATTEMPTED computation inside it (a `try/except`
  that prints the exception type), run on each discriminating configuration;
  `inspect.signature` bounds the reach of a hole.
- **L-026** "Is this degraded condition deliberate?": run the test at the
  healthy setting and re-evaluate its OWN assertion. The unit is the SOLVE,
  not the test (the per-test call count predicts a second solve). Never cost
  a budget raise by the budget ratio.
- **L-016** A test asserting `== <the stored tag>` pins the tag and is
  evidence the property is UNTESTED. When a flag's docstring and its registry
  entries disagree, measure which is wrong.

Sharpening the definition §3 (an inherited datum is a claim):

- **L-011** A docstring that DELEGATES ("the sweep handles it"): grep the
  mechanism NOUN the claim invents, never the symbol. A sibling
  implementation that REFUSES the same input is a free oracle, and a strict
  `xfail` row naming the gap is the tree's own admission.
- **L-021** A brief's arrow `A --f--> B` is a claim that `f` is a
  MATERIALIZED object: grep its accessor's return type before mapping sites.
  An earlier phase's survey of the same campaign expires like any other.
- **L-023** "N spellings of one concept": read the CALLEE's `return` before
  mapping the callers, since the N sites are usually N callers re-deriving a
  fact the callee computed and dropped. Grep the DEFAULT of a boolean claim
  field. "Raise or warn?": grep the pytest configuration for
  `filterwarnings` first. (§L-020: check that the spellings share a domain
  and codomain before endorsing a unification.)
- **L-020** When re-auditing a negative claim, read the NEW module that an
  intervening commit added and ask the claim against it.
- **L-027** Reconciling a survey: `git diff --stat <baseline>..HEAD` on the
  PRIMITIVE the survey is about, and read its `__post_init__`; a new required
  field is an obligation no claim-by-claim check sees. A handed count is
  printed under several conventions side by side, and its qualifier ("all in
  `tests/`") is checked as a separate claim.
- **L-035** A claim confirmed at its cited site is unconfirmed until the
  site's CONSUMERS are counted; the verdict that carries both facts is
  "CONFIRMED (substance) / REFUTED (site)".
- **L-036** Zero commits since a plan's date proves staleness INHERITED, not
  drift. A plan that prescribes what a ratified counter-design replaced is a
  "conflict to adjudicate", neither NOT-LANDED nor REFUTED.
- **L-034b** A campaign that ships its own gate suite: run it with `-rx`;
  the `xfail` rows are the premise oracle.
- **L-012** A carve already underway: diff the uncommitted hunks against the
  brief's items (one may be mid-fix), reconcile each finding by a runtime
  probe against the final tree and not by reading the diff, and lead with
  the items the carve did not reach.

Sharpening the definition §4 and `retirement-audit` (blast radius):

- **L-017** The frozen baselines a retirement can move are found by
  `find tests -name '*.npz' -o -name '*.npy'`, never by a symbol grep
  (`retirement-audit` D.14 carries the independence question).
- **L-028** A units change is often an ORDER change, and may put a ZERO in
  the new range: ask which denominator reads `edges[0]` or `edges[-1]`, and
  whether any test asserts the producer's returned VALUES at all.
- **L-038** Print `__all__` against the AST public set: the diff surfaced
  production-consumed names missing from it.
- **L-040** A guard-message pin census maps each `match=` to the OPERATOR it
  exercises. A step named by a code is checked for homonym codes before the
  census. A chartered assertion with 0 spellings in the tree means the step
  OWES the gate: say so, because the zero reads as "nothing to do". A role
  whose explicit target is the CLASS is never flagged dead after its method
  dies.
- **L-031** Citation drift from a read window is uniform, so sampling three
  citations detects it.
