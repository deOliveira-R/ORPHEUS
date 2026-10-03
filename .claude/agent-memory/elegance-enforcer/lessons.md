# Elegance Enforcer — Lessons (hot digest)

Review-PROCESS corrections only: "what did I mis-judge or catch, and what sharpened
the verdict?" Read this every dispatch.

**Scoping — check these three owners before adding anything here.**
- The pattern/anti-pattern CATALOG is the preloaded `coding-elegance` skill. NEVER
  restate it; cite it (`Pattern N`, `anti-#N`). Several lessons that used to live here
  were PROMOTED INTO the skill — anti-#20 ⊃ L-001, Pattern-4 `replace()` ⊃ L-008,
  Pattern-6 self-conceding-trim ⊃ L-007, anti-#19 ⊃ L-012. What survives below is only
  the *detection* half: how the smell LOOKS in a diff, and how I verify it.
- The institutional SMELLS (twin-delivery plumbing, role grid, fuller-view oracle,
  tells-to-grep) are AGENT.md §"Institutional knowledge".
- The three-leg VIOLATION standard (name the future edit / coextensive-today→NIT /
  verify the LIVE tree) is AGENT.md §"The VIOLATION standard". It governs every
  verdict by definition; each lesson below is one face of it.

**War stories, `file:line` forensics, per-review inventories → `lessons_archive.md`**
— the byte-identical predecessor of this file; `→ archive L-NNN` points at its headings.
L-NNN ids are load-bearing: sibling memo files cite them.

## Standing review order (every dispatch)

0. **⭐⭐ ADVERSARIAL PHASE FIRST — see L-020.** Never open with a balanced
   survey. Ask "how would I BREAK this" (silent wrongness, weighted to defects
   that fail in the REASSURING direction) and "how would I make this 100×
   better" (reframe the JOB). Balance, and every "well-factored / do not touch",
   is Phase 2 — written as a WITHDRAWN ATTACK with the reason you expected a
   defect.
1. **Scope from the LIVE tree, never the brief.** `git status --short` + `git diff
   --stat HEAD` FRESH. Untracked (`??`) files never appear in `git diff` and are the
   easiest coverage to miss — then grep the tested SYMBOLS across `tests/`. (L-011)
2. **A surprising NULL is a tooling bug until re-verified with a differently-spelled
   grep.** This self-catch fired three times and each time was about to become a false
   MUST-FIX: zsh glob-expanded an unquoted `--include=*.rst`; a `\|` was literal under
   `grep -E`. Never flag an absence off one pattern.
3. **Filter `_build` on any tree-wide sweep.** Gitignored
   `docs/_build/html/_sources/**` carries OLD anchor definitions and OLD page paths —
   an apparent twin that is build cruft. `| grep -v _build`, or `test -f` the source.
4. **Grade with the three legs.** If you cannot name the future edit that makes the two
   things diverge, it is a CONCERN or NIT — write it as one.

## A. Verify before you flag (leg-3 machinery)

### L-012 — PROVE a `# type: ignore` is dead; the error-count ratchet is BLIND to it
`reportUnnecessaryTypeIgnoreComment` defaults OFF, so a dead ignore survives every
error-count ratchet — and two carves CLAIMED retirements the live tree had not made.
Don't eyeball: write a throwaway config under `scratch/` mirroring
`pyrightconfig.json` + that rule as `"error"` + `include` the one file, pass it with
`-p`, grep the lines (the write-scope hook refuses a config at the repository root). Read-only — never edit the file
under review. Run it on any file that GAINS **or MOVES** an ignore: a relocated ignore
can have been dead before the move, and a precedent-setting carve makes the dead ignore
the template every sibling copies. → archive L-012

### L-017 — "This bound rejects that carrier" is a CONSUMER-site claim — CLI-pyright it
I built a near-certain SHOULD-FIX by reading a `V bound=Vector` class-generic beside its
ndarray-driving L0 tests. `npx pyright` on the consumer REFUTED it: a class-level TypeVar
solves PER-CALL, and an unpinned `V` swallows the mismatch — zero errors at the
instantiation. Never infer a rejection from the TypeVar declaration; the definition file
is almost always clean. Run pyright at the CONSUMER and attribute every diagnostic to a
line. The defensible verdict on a typing carve is earned by production-file pyright-0 +
a suppression grep, never by the docstring's self-description. → archive L-017

## B. Grading the finding

### L-003 — Parallel predicates with different return types are NOT a unify trigger
They speak to different consumers at different layers. Twice asked "these two gates look parallel — co-locate them?" and NO was right. A `bool`
physics-validity query and a `Compatibility(ok, reason)` dispatcher query are parallel in
spirit but sit on different axes at different layers. If unifying forces a boolean flag
(anti-#3), a cross-family dependency, or a `Compatibility`-vs-`bool` coercion, it is a
CONCEPT MERGE — the opposite of a single-source win. Keep them apart and say why. Same
logic: a CONSTRAINT that defines an operation's identity (convex vs affine combination)
is not a flag parameter. → archive L-003

## C. Blast radius — outside the diff, and inside the edited file

### L-016b — A LANDED capability inverts the blast radius onto the stale DEFERRAL CONTRACT
The mirror of L-004: when a carve flips deferred→implemented, every docstring naming the
case "raises / deferred" is present-tense-FALSE. The half-cleanup signature — which
recurred even when the author explicitly applied this lesson up front — is: the
human-facing **rst ledger** gets rewritten, the machine-facing contracts do not. The four
sites that keep surviving: the `@runtime_checkable` **Protocol stub**, the **base default**
docstring the next implementer inherits, the **sibling operator CLASS** docstring, and
**public operators in files the diff never touched**. Grep the deferred case's name AND
its prose forms across the whole package, then discriminate BY ARM/SET: a matvec-transpose
landing does NOT un-defer the transpose-SOLVE; only the adjoint row flips while α/transient
rows stay genuine future seams. In a campaign-CLOSING phase it is MUST-FIX — "closing #NNN"
is internally inconsistent with a row still tagged "future seam". Milder siblings: a
blocker that CLEARED while the deferral stayed true (reason-staled ⟹ SHOULD-FIX); an
in-diff comment describing a plan the same diff's later commit executed (L-015 in-diff).
→ archive L-016 (2nd heading, + its three sharpenings)

### L-016a — Certifying the ARRIVAL of parked scaffold (the inverse of a gather/split)
Enumerate every parked atomic fact from `git show HEAD:`, bucket each {present |
ADAPTED-as-declared | LOSSY}, and OPEN the claimed new home for any declared drop — one
"loss" was a Pattern-2/7 WIN (a solver-specific magic number leaving a cross-cutting
doctrine page for its native chapter). Char-diff the declared-verbatim blocks. **The
killer one-home check is STRUCTURAL, not the anchor grep**: a rival definition rarely
reuses the moved label, so grep the concept's SHAPE (the list-table row form, the
"L0 …, L1 …" enumeration). And the surviving twin hides in the SAME file as a repointed
site, outside the diff — read the whole owning file hunting a second definitional verb
("is defined"). Doctrine-page-twins-a-skill is CONCERN-not-VIOLATION given reciprocal
pointers + an explicit ownership rule + a functional partition; the residual habitat is
that the ownership rule is unenforced prose. → archive L-016 (1st heading)

### L-018 — A class SPLIT strands an IN-DIFF docstring contradiction — turn the lens INWARD
After the outward `git grep`, re-read the EDITED file's own module/class docstrings top to
bottom. The split INSERTED a correct "Two discrimination axes" section spelling
`SolutionBase.is_eigenvalue` and left the intro 20 lines above still spelling
`Solution.is_eigenvalue` — one docstring, two spellings of one fact. The recurring
LLM-refactor shape: ADD corrected prose next to the change, never reconcile the older prose
the change just falsified. In-diff ⟹ SHOULD-FIX now, unambiguously in scope. Leg-3
restraint that kept me honest: with no build available I cited the PROSE drift (certain
from reading) as the finding and flagged the `:meth:`Leaf.moved_method`` xref break as
"verify on next build" — never upgrade an unverified xref break into the finding itself.
→ archive L-018

## D. Elegance calls (where the skill's catalog meets a judgment)

### L-007 — Grep for PRODUCTION consumers before crediting any new predicate/Protocol/ABC
Catalogued as Pattern-6's self-conceding-docstring trim. Detection: zero production
consumers + a docstring conceding production won't use it ⟹ TRIM on sight, and retarget
its tests to the pre-existing lower-level property it was minted to pin. KEEP only if a
real production consumer lands in the SAME commit. Contrast I nearly over-called: a
dunder-EMPTY role-marker mixin whose second consumer is IN THE SAME DIFF satisfies
rule-of-two at landing — and a long docstring on an empty body is ELEGANT there, because
the content is "absence of a gate," which cannot self-document via code. → archive L-007

### L-008 — Prefer the emergent invariant-gate over a hand-written dunder
Catalogued in Pattern 4 (route same-type-producing ops through `replace()`). The review
move: asked "should this type forbid operation X?", check whether a construction invariant
already exists that `replace()` will RE-RUN. If yes, an explicit `__add__` override is a
Pattern-2 duplicate of the law — flag it. If NO invariant exists (flux+flux is
type-coherent but void), the explicit gate is CORRECT. Caveat worth a clarifying comment:
`@runtime_checkable` Protocol conformance is method-PRESENCE only, so a gated-arithmetic
leaf passes `isinstance` and then raises inside a generic accumulator loop.
→ archive L-008

### L-010 — On a library/tooling review, DIFF THE GUARD CLAUSES of sibling methods
A guard present in N−1 siblings and absent in 1 is the latent bug (symmetry-in-code). The
Nexus case: one overlay method omitted the `if node_id in self._g` stale-node guard, and
the crash is SILENT — `degree()` on a missing node returns a view, not an int, so it
surfaces only as `TypeError: not JSON serializable` at the MCP boundary. Also flag the
self-disproving docstring (an "exactly one metric family" invariant refuted by the same
commit's heterogeneous-union merge): "it's just a docstring" makes it MORE insidious, not
less — it is the contract the next agent reads. → archive L-010

### L-019 — A builder↔display WELD can still harbor an unwelded THIRD inline spelling
In an algebra-of-record (`derivations/`) module, do NOT credit the weld from its docstring.
Grep EACH theorem body for whether it CALLS the builder or RE-SPELLS the rule inline. A
`_display_matches_builder` helper welds only the pair it names; a downstream theorem's
local `t1()/t2()` byte-identical to the canonical builders is an unwelded twin — the day
the builder is refined, the theorem proves a STALE rule and still passes. Discriminator:
a rule a builder already provides = twin (SHOULD-FIX; fix by calling the builder) vs a
deliberate counterexample = legit. NOT a twin: production numpy vs the SymPy builder —
that is the intended derivation/impl structural independence. → archive L-019

## E. Doc-carve certification (my largest recurring workload)

### L-013 — A doc gather/split/move is a RETIREMENT review; mechanize it, never eyeball
Sixteen chapter-carve reviews collapse to seven rules. Instances (all 2026-07,
COMMIT-READY or 1 MUST-FIX): promotion / multi-span-differential-depth / chapter-MINT /
DEMOTION / wrapper-dissolution / H1-dissolved-3-ways / fresh-authored chapter / router /
NO-CHANGE adjudication / notation-harmonization / page-level `git mv`. → archive L-013 +
Sharpenings 1–16.

1. **The certification kit.** Per-span `difflib` char-diff against `git show HEAD:<page>`;
   full-file `@@`-hunk enumeration (a per-span diff is BLIND to an out-of-span edit);
   underline/promotion census by type with LENGTHS preserved; header tree extracted and
   asserted for ZERO level jumps (classify PER SPAN — a union map is ambiguous when `~`→`=`
   in one span and `~`→`-` in another); label single-homing; rump checked for orphans.
   Account for EVERY hunk against {declared transform | traveled | merge | seam-tidy}.
   Declared in-span fixes are the expected residual, and each must repair a GENUINE
   falsehood — grep the referenced content's ACTUAL home; a fix that maps to no real home
   is editorializing.
2. **FRESH vs TRAVELED is the master verdict discriminator.** Settle it by char-diffing
   the paragraph against HEAD. Fresh-edited prose carrying a stale/incorrect spelling ⟹
   MUST-FIX. A byte-exact traveled span carrying the same spelling ⟹ CONCERN, deferrable
   to the declared harmonization stage. Hold ONLY fresh content to the current standard —
   UNLESS the author touched it, which pulls it back into scope.
3. **The live finding is ALMOST ALWAYS the external blast radius (L-004), and the author's
   blind spots are `tests/` and the FOUNDATION docs.** Causation oracle: `git show
   HEAD:<old page> | grep -c <label>` NONZERO **and** worktree `grep -c` ZERO ⟹ THIS stage
   staled it ⟹ MUST-FIX; both nonzero ⟹ pre-existing ⟹ SHOULD-FIX-for-completeness. Run it
   PER LABEL, because a compound reference can cite one moved and one retained label — and
   that partition IS the MUST-vs-SHOULD line (fully-transferred ownership = MUST; a
   half-stale page attribution whose sibling label still anchors the prose = SHOULD). A
   PARTIAL test sweep is more insidious than none (the diff SHOWS tests being repointed);
   grep the bare moved path tree-wide and demand ZERO. A page-level `git mv` upgrades
   severity: the `:doc:` is BROKEN, not merely qualifier-stale. `@pytest.mark.verifies` is
   page-agnostic — never flag it; only prose/docstring page attributions.
4. **A fresh-authored page (router, thin chapter, corpus root) has no char-identity
   surface — the whole text is the claim-truth surface, and the recurring catch is an
   OVERCLAIM.** Spot-check every load-bearing claim against the LIVE code; a verification
   section is a table of grep-checkable assertions — check them ALL (test names,
   tolerances, xfails, guards). Grep strong identity words ("byte-for-byte", "identical",
   "the same … verbatim"): the proof it is an overclaim is that the chapter states the
   weaker TRUE version of the same claim elsewhere in itself. For a ROUTER, the product is
   its forward links — `-W` proves a `:doc:` RESOLVES but not that the target CONTAINS the
   promised derivation; extract the target's header tree. And if a SIBLING page makes the
   same claim, the gap is corpus-INHERITED ⟹ CONCERN + issue, not a fresh MUST-FIX.
5. **The project spelling axis.** Honest algebra is `A = L + C − S − B` (A = the FULL
   within-group operator); the sweep is `(L+C)^{-1}`, the inner kernel of the full inverse,
   NEVER the full inverse itself — so `A = L+C` and `A^{-1}`-for-the-sweep are the stale
   spellings to grep for, inline `:math:` included, not just `.. math::` blocks. A LOCAL
   `A` explicitly bound at point of use is the strongest honest form, not a defect; the fix
   for a cross-layer notation seam is the IN-CODE bridge that spells both bindings, NOT
   forcing one letter. Fission is cross-group — `−F` in a within-group operator is a
   physics error (Cardinal Rule 1), not a notation NIT.
6. **On a notation-harmonization pass, grep the PROSE forms separately from the formula
   forms.** The author's scan matches `A = …` / `A⁻¹` and misses "the operator :math:`A`"
   in the summary bullet and the section intro that PARAPHRASE the renamed formula —
   a pass-introduced internal inconsistency. And respect the math/code seam: rename the
   math symbol, but a double-backticked ``A`` names a CODE construct — renaming it to a
   symbol the code doesn't have drifts the doc from the code (anti-#20).
7. **Adjudicating a NO-CHANGE / anti-padding ruling** (a distinct review shape): run a
   genuine ADVERSARIAL pass first — enumerate the candidate additions and test each — then
   concur. The decisive confirmation is that the proposed addition would create a TWIN, not
   merely be padding. Also: a router routes on its natural axis (concept map vs task table
   vs symptom table); forcing a second axis onto it ADDS a concept without collapsing one.
   A corpus root ORIENTS-AND-DELEGATES; smuggling a part-index's routing table upward is
   the violation to check for.

### L-014 — Certifying a code-prose REBALANCE (docstrings trimmed to theory-book pointers)
Recipe in `doc_prose_rebalance_certification.md`; forensics in the archive. The spine:
**behavior-invariance = token-invariance (drop COMMENT+STRING+layout) AND a
docstring-stripped `ast.dump()` compare** — the AST leg is strictly stronger, since it
still SEES non-docstring literals, so a change hidden in an error-message or dispatch
string cannot slip past the STRING filter. (Py3.14 gotcha: f-string literals tokenize as
`FSTRING_*`, not `STRING`, so they ride in the code view.) Pointer honesty = resolve THEN
content-check the landing site (leg-3). The MUST-FIX class is **contract
self-sufficiency**: a caller must work without leaving the file (mutation semantics,
producer conventions, raise conditions). Calibration: a machinery/driver file yields a
~14× smaller cut than a teaching file — CORRECT, not under-delivery; hunt the `#`-comment
retirement TOMBSTONES there, and grep each cut tombstone's claimed destination to confirm
the constraint landed. → archive L-014

### L-021 — A registry key added ONLY to route to a REFUSAL corrupts every consumer that reads the registry as an INVENTORY
Detection, and it is grep-then-RUN. When a diff adds table/registry keys whose declared
purpose is "reach the handler and be refused there" (a decoy key giving a better error than
"no entry"), the fix is local and the damage is not: **find every site that enumerates the
registry's keys and run it.** The recurring victim is the sibling error path — a
`NotImplementedError` that prints `sorted(REGISTRY)` as *"catalogued today"* now advertises
entries that unconditionally raise, so the message meant to be the map of what EXISTS is the
one place still selling the retired spelling as live. `[M]` 2026-09-02, #432: 3 decoy
`Sphere/SO2_a` keys made `SPHERE.quotient(Cn(4))` list 9 entries of which 3 are unobtainable.
Grading: this is NOT a coextensive-today NIT — *catalogued* and *obtainable* disagree the
moment the key lands, so leg 2 of the VIOLATION standard is met without any future edit.
Structural fix to demand: split the two roles (a `_ALIASES` table consulted before the
registry), never a filter in the message — a filter is the same fact spelled a third time.
⭐ Generalises §10's shape (a metric invalidated by its own campaign's success) from a NUMBER
to an ENUMERATED SET: ask of any registry-derived listing, *"after this change, is every
member still deliverable?"*

### L-022 — A refusal placed in the DERIVATION instead of on the TYPE pays in import edges and re-derived arguments
The tier question ("where does this new invariant live?") has a cheap tell: count what the
chosen site had to IMPORT and RECONSTRUCT to state it. `[M]` #432 put "the `by` group must be
the stabiliser" inside a catalogue derivation, which then needed (a) a new function-scope
`manifold -> symmetry` runtime import in a module whose docstring is a measured argument about
that exact cycle, (b) a `letter = "xyz"[axis]` round-trip because the accessor it had
(`rotation_axis: int`) had thrown the letter away, and (c) decoy registry keys to be reachable
at all (L-021). All three vanish if the fact lives as a property on the GROUP and the check on
`__post_init__`. ⭐ The confirming signal is usually already in the file: here `Quotient.__post_init__`'s
own docstring says *"a mis-specified entry is refused **where it is written**, not where it is
read"* — the doctrine was present and the new invariant did not follow it. So: before grading a
refusal's tier, grep the target type's existing `__post_init__` for a doctrine sentence, and test
`dataclasses.replace(entry, <field>=<illegal>)` — Pattern 4 promises `replace()` re-runs the
invariant, and it only keeps that promise for invariants that are actually IN `__post_init__`.

### L-023 — A field excluded from `__eq__` + a memo keyed on the owner = a CALL-ORDER-dependent answer
Detection, and it is a two-step grep. When a diff adds a field whose stated job is to let a
consumer DISCRIMINATE ("the codomain a reader checks to learn which it was handed", a
`kind`, a `units`, a `provenance`), ask (a) can `__eq__`/`__hash__` tell two owners apart on
it, and (b) is there a `functools.cache`/dict keyed on the OWNER? Both together and the memo
serves one owner's answer for the other's — and the corruption flows from the forged object
INTO the honest one, so it is not merely "the liar lies". `[M]` 2026-09-03 #434 R4:
`barycentre(entry).codomain` reads `S^2` or `D^3` purely by which entry was asked first.
⚠ Grade this a VIOLATION without a future-edit hunt: leg 2 is met the moment the field
lands, exactly as in L-021 (*catalogued* vs *obtainable*). The tell in a diff is a
`field(compare=False)` whose justifying comment says "derived from X" while a SIBLING field
derived from the same X is compared — that reason proves too much; the honest reason for
excluding the siblings is usually "a function has no value equality", which does not
transfer to a value-typed field. → topic file `symmetry_realization_carve_rulings.md`

### L-025 — On a DESIGN review of a type-tightening, count the producers of the refused state at RUNTIME; the guard's own docstring is the producer index
Pattern 6's guardrail says "grep for a production path that legitimately produces it"; on a
W5 design (no code yet, so no pyright) the grep is uninformative when every constructor site
spells the field the same way. `[M]` 2026-09-27, posing-sequence attack: the section ruled
"`SigT` becomes a property, the guard retires, conservation by construction"; all 8 production
and 32 test `Mixture(` sites spelled `SigT=` by keyword, and the guard `assert_balanced` was
deliberately OUTSIDE `__post_init__` with a docstring naming three legitimate imbalanced
producers. A `-p` pytest plugin wrapping `__post_init__` over the files that spell the
constructor (23 files, 6 min) read **167 of 584 constructions imbalanced**, one class in
PRODUCTION derivations (a benchmark registry encoding), and the docstring's own account of that
encoding was stale. Detection order: (1) read the guard's docstring for the reason it is not
in `__post_init__`, since that reason is the producer list; (2) wrap the constructor's invariant
in-process and count over the constructor-spelling files, attributing by the frame ABOVE any
shared factory; (3) name the re-spelling of each producer class in the finding, because the
design's DIRECTION (one definition, Pattern 7) is usually right and the census is the only
missing leg. Grade VIOLATION when a production producer exists (anti-#18 leg (ii)).

### L-024 — Tests that RE-COMPOSE the SUT's formula from its fields pin the tests' formula, not the module's
Not covered by `retirement-audit` D.14 (demotion by a *retirement*) nor X4's tells
(`allclose(a, b)`): here nothing was retired — the SUT's composite property is too slow
symbolically at one size, so each test file re-spells `R A⁻¹ E + f I` pointwise from the
SUT's fields, and the published property ends with ZERO pins wherever no in-module
proof reads it. Detection: for every composite property of an algebra-of-record type,
grep the tests for the property NAME (not the class); if they read only the fields,
mutate the property in-process via a `-p` pytest plugin (+ a positive control on a
covered arm). `[M]` 2026-09-23 face transmission: `CellResponse.transmission` negated for
LD d≥2 → 0 of 83 red; control → 4 red. Remedy is one composition site with a pointwise
entry (`transmission_at`), never a third spelling. Also: check the slowness claim that
justified the re-spelling — it was stale (0.08 s vs "does not finish in 30 s").

### L-027 — On a PLACEMENT review (where does a datum live: basis / measure / coupling), prototype every placement as the SAME numpy chain and census the READERS of what the object builds
Not covered by L-025 (producers of a refused state) nor the brief's own census ask. `[M]`
2026-09-28, the PG-weight reversal: (1) an AST census of `.basis_space/.test_space/
.measure_space` readers showed 0 of 7 weighted frames have a downstream reader — every
"space identity / caching / `.H`" attack on the placement was then either withdrawn or
re-graded pre-existing (F5) in one step, instead of argued; (2) prototyping A/B/C on one
fixture told bit-identity (C ≡ A, 12 of 12) from float order (B, 1e-16 on 7 of 12), which
is the whole answer to "do the 0-ULP pins survive"; (3) a `-p` plugin wrapping the carrier's
constructor over the consuming suites (177 tests, 9 s) measured the refused state's runtime
population (0 of 106 negative), the L-025 move applied to a POSITIVITY law; (4) a spy on
`Basis.evaluate` counted one table tabulated 4×/7× per verb — the cheapest "does this
placement re-do work" instrument, and it found a Pattern-2 twin (`membership` beside three
frames) the read had missed. Also: a wrapper `Basis` that raises on half its surface makes
the FACE's `is_adjointable` role default (`True`) a representable lie — test `.H.apply`,
never read the flag. And read `git log` for the ruling's history before grading a reversal:
this one had already flipped once on an unmeasured ground (F8), so the finding is "land it
with a witness", not "land it".

### L-028 — An AST Name census is BLIND to a quoted forward-reference annotation; and on a RENAME set, census the DESTINATIONS first
Sharpens standing order 2 (a surprising null is a tooling bug) by naming the spelling that
hides: `def f(m: "TransportMethod[OpT]")` is a `Constant` string with brackets, so a
`Name`/`Attribute` walk and an identifier-only string filter both report 0 consumers, and I
had a Pattern-6 "zero production consumers, trim" finding drafted before a `def` grep showed
the one consumer. `[M]` 2026-09-28 posing-sequence W5. Before any zero-consumer verdict,
grep `["']Name\b` and `Name\[` as well. Second half, not in `retirement-audit` F.23 (which
covers a rename's RESIDUE): a design that renames N symbols owes a census of the N
DESTINATIONS for existing occupants holding a DIFFERENT quantity; `[M]` 3 of 24 were taken
(`balance_defect` = a relative norm, `ordinate_idx` = the global index the cylinder arm does
not pass, `provenance` = two citation-record classes), and each was the review's top finding.
Also from that review: prototype a design's "bit-identical" claim on the CORNER fixture (a
straddling grid), not the nested one; it is where the retyped verb and its would-be twin
(`(MR)⁻¹M` vs the row-sum ratio) separate (0.36) while agreeing everywhere else.

### L-029 — On a design that "collapses an enum into a derived pattern", count the enum's production DISCRIMINATION sites first; on a control battery that permutes the design's own data, feed it an ALIEN from the other fixture
Not covered by the skill's "a repeated conditional is a missing type" (which presumes a
conditional exists) nor by X1 in general. `[M]` 2026-09-28, the coproduct W5: my AST dispatch
census printed 0 production sites with a FAILED control (the control I chose, `a if a is b else
Enum.X`, is a `Compare` with no enum operand) — a grep of every line naming a member, with the
join line as control, then hand classification, found the same 0 branches, so the enum was
test-pinned metadata (44 of 44 `isinstance` marker reads in `tests/`), and the finding was L-007's
trim plus a present-tense-false "dispatches on" comment, NOT a missing type; grading it as the
memo framed it would have credited the design for retiring a dispatch that did not exist.
Second half: 48 mis-ORDER controls all reddened and none could see that the new arrows admit a
composite from another mesh (production's `admit_composite` refuses it); the alien was the other
fixture the battery already held. Third: reproduce a concept-count claim ("11 → 8") under a
STATED definition and add back what the memo's own text says it KEEPS (28 → 31 as prototyped).

### L-026 — On a "every X is the <induced arrow> of one map" unification, check the ROLE the arrow plays per instance
Not covered by the type-vs-property test (`coding-standards`), which a real concept passes
while the role assignment is still wrong. `[M]` 2026-09-27 (`Pullback(φ)` W5): the pullback
is the RESTRICTION along an injective map (trace, system member) and the unnormalised
EXTENSION along a surjective one (axis projection: `R = (π*)^H`, `E = π*/Σw`); where the Gram
factor sits is a physics convention (ERR-051) no map can derive. Detection: per instance,
ask whether φ is injective or surjective and which induced arrow the tree calls the
retraction; grep the theory corpus for a sentence already refuting the identification
(`spaces.rst` "The pullback is not the section"). Also probe a point-map type spelled by
node POSITION (`lambda nodes: nodes[idx]`) on a reordered array: it answered wrongly.

### L-030 — On a GRAPH-structured design (an ordering, a cut, an SCC), construct the input whose declared effect the STRUCTURE absorbs, and run it to the answer; and check whether a declared-order constructor DERIVES its cut or makes the call site compute it
Not covered by L-025 (producers of a refused state), L-029 (an alien from the other fixture) or the
mutation battery, which all vary the design's INPUTS; here the state is a legal input whose effect
the structure swallows. `[M]` 2026-09-29, the `Ordering` W5: every fixture cut broke its SCC, so 13
gates and 4 mutation arms were blind to a cut whose ends stay strongly connected through a second
cycle — accepted, lagged AND kept inside the piece's own block, "converged" in 8 passes to 7.8e-2
(control 1.7e-15). Build the degenerate by hand (a 3-block second cycle), solve to the answer, and
read the piece's own block for the cut entry. Second half: production's Gauss–Seidel DERIVES its
lagged rows from the declared order (`lower_inflow_rows`); the prototype REFUSED an order with
upstream edges and its task file computed the cut with three lines of index arithmetic above the
"call site" marker — grep the prototype's task file for arithmetic between the fixture and the
marker, because the memo counts statements only below it. Third: each "made unspellable" claim in a
prototype memo is a prose claim (X3); spy the verb it names (`Split.M` was called once on the very
path the memo called unassembled).

### L-031 — A diff that adds a ledger TOKEN owes the tree-wide ledger gate a run; and feed every ALTERNATE constructor the inputs the primary's parse refuses
Not covered by the definition's §4 instruments (they check a gate's teeth, not whether a new token satisfies an
existing global gate) nor by Pattern 4 (which states parse-don't-validate, not how to find the bypass). `[M]`
2026-09-29, P1 step 2: the new `SCOPE-BOUNDARY[guard]` put `ruling:`/`revisit:` outside
`test_elegance_debt_is_tagged.py`'s 7-line window. It was red on the working tree, and the step's own gate
files did not include the ledger. So grep the diff for `SCOPE-BOUNDARY|ELEGANCE-DEBT|verifies\(|catches\(`,
and run the gate that reads each token. Second: a vocabulary classmethod (`from_thicknesses`) that pre-coerces
with `float()` before delegating admitted `"0.5"` and `True`, which `__post_init__` refuses. Probe each
`from_*` with the primary's refusal fixtures, in-process.

### L-032 — A twin kept "for bit identity" is graded by replaying the regression CAPTURE under the collapsed spelling, and by asking which spelling is correctly rounded
Not covered by leg 2's coextensiveness check (the two spellings DISAGREE, and the disagreement is the stated reason to keep
both) nor by L-028 (prototype on the corner fixture). `[M]` 2026-09-29, #405 P1 3a: `interval_measure` (Python `b**2` = libm
`pow`) was kept beside `compute_volumes_1d` (numpy square) "for ERR-020 bit identity". Against a `Fraction` oracle, the kept
spelling was the INEXACT one (22 of 20 000 against 0), and replaying the pre-carve capture's 414 curvilinear equal-volume
intervals under the array spelling moved 0 (control: the re-association moved 6). The identity protected only a gate against a
function retiring next step. So: find the capture the identity claims to protect, replay it with a positive control, and
test both spellings against an exact oracle before crediting "kept for bit identity".

### L-033 — An admission keyed on a family TAG: read the tag's definition first, then feed the guard every PARAMETER of the family, using a sibling method's realizer as the oracle of what is ill-posed
Not covered by L-031 (an alternate constructor of ONE type against its primary's refusals) nor the skill's anti-#4 (which
names the smell, not the probe). `[M]` 2026-09-30, ERR-094 review: I had "CP's `kind == "reflective"` guard re-admits a
partial reflector" drafted; `ReflectiveBoundary.kind` returns `"partial"` at α ≠ 1 (a hidden fourth `albedo == 1.0`), so it
was refuted by one read. The real hole was the parameter the bug report did NOT name: `ReflectiveBoundary("y")` on an
x-face was admitted by CP (k bit-equal to the x-mirror) while SN's realizer refuses it as "not a boundary law at all".
So: (1) read the tag property's body before grading what a string guard admits; (2) enumerate the family's parameters
(axis, sign, amplitude) and run each through the guard AND through a sibling method that realizes the law, whose
refusals are the free catalogue of ill-posed members; (3) the remedy is value equality with the law the method computes,
not a structural predicate (`law_permutes_ordinates` also answers True for the wrong-axis mirror).

### L-034 — When a design RE-HOMES a type to a new consumer, feed the new consumer the values the type stores today; guards written for the old consumer's blindness strip what the new one needs
Not covered by L-031 (alternate constructors against one type's refusals) nor the skill's Pattern 4 (which asks what a type
admits, not whether its stored data suffices for a different reader). `[M]` 2026-10-01, orbifold prototypes W5: decks were
to move from the face law (read by the realizer, which sees only the linear part) to the geometry (which must generate the
deck group). `SelfPairedDeck` refuses a mirror's offset and `PairedDeck` stores a unit wrap, both correctly, as Mode-12
closures for the realizer. Fed to the geometry's consumer, `close_group` of a both-faces-mirrored slab's stored decks
returned order 2 (Z2) for the true D-infinity, and the two decks were equal values. So: name the new consumer's verb,
run it on the stored values from a shipped fixture, and expect the data model to invert (store the located object; derive
the old consumer's view from it).

### L-035 — A design claiming "the key covers every class's SCHEMA" is graded by diffing, at runtime, each class's emitted part names against its declared fields; a hand-enumerated part list is a second schema the claim and its gate are both blind to
Not covered by L-028 (a census blind to quoted annotations) nor Pattern 2 as stated (it names the twin, not the probe that
finds it). `[M]` 2026-10-02, #405 P1 step 5: the encoder's schema tag is the part NAMES, and `Axis` overrode
`content_parts` to drop `generator` by omission while `Mesh1D` used `field(compare=False)`. So a field added to any of 4 of 22
classes entered neither the digest nor the tag, and the S5.3 "a part added later reds" gate read the same override. So:
(1) list `dataclasses.fields` against `content_parts()` names for every subclass (walk `__subclasses__`); (2) any exclusion
spelled by omission rather than at the field declaration is a VIOLATION, remedied with `compare=False`; (3) ask which
ladder-order dependence the override was hiding (here: ContentIdentity-before-dataclass).

### L-036 — A verdict built from several PAIRWISE consistency checks over enclosures of ONE exact quantity is graded by Helly (1-D): build the state where the UNCHECKED pair is disjoint; and an "exact" value evaluated by `evalf` is probed with a cancellation-to-zero expression
Not covered by Pattern 2 as stated (the twin here is a missing LAW, not a duplicated body) nor by L-030 (graph absorption).
`[M]` 2026-10-03, #405 P2 step 5: `ReferenceCertificate.state` checked claim∩anchor, refinement members together, and
claim∩common part; anchor∩refinement was never asked, so a 0.50 anchor and refinement members at [9.2, 10] read `Valid`.
A family of intervals shares a point iff every pair meets, so the law is ONE intersection over all enclosures per
observable. Same review: `Exact` of `cos(π/7)−cos(2π/7)+cos(3π/7)−1/2` (exactly 0, not auto-simplified) enclosed
−1.4e−191 ± 1.4e−250, excluding 0; `evalf(n, strict=True)` raises `PrecisionExhausted` instead. So: list the pairs a
consistency verdict checks against all pairs of its inputs; feed an "exact" evaluator an identity that cancels to zero.
