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
