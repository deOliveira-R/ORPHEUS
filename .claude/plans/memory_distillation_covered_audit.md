# qa audit of the COVERED claim class, memory distillation 2026-11-16

Subject: every row an owner marked COVERED in `scratch/memory_distillation/<owner>/dispositions.md`.
Worktree `/Users/rodrigo/git/nuclear/ORPHEUS_mem`, branch `chore/memory-distillation`, HEAD at start `bcd4478b`.
Carriers read in the worktree (`git diff 0d68fc28 HEAD -- .claude/agents .claude/skills .claude/rules` touches only
`code-search.md`, `workflows.md`, the `instrument-doctrine` skill's ugrep sentence and two workflow-stage agents; no cited clause moved).

Method, per row: the retired body was read at its pre-retirement commit (for a lesson retired on 2026-09-22 the parent of that
retirement, NOT `1039683f`, where those lessons are already one-line stubs), its distinct checks listed (last paragraph
included), each check located at a carrier SENTENCE the owner loads (front matter `skills:` and `omitClaudeMd` read per owner),
and any miss searched in the owner's worktree memory directory. Census: a Python `re` script over the named files, with a
positive control per owner (a token known present).

Verdicts: CARRIED (every distinct check is in a loaded carrier or kept in the owner's new memory); PARTIAL (one check is in no
loaded carrier and in no file of the owner's new memory; the check is quoted and a restoration text given); FALSE.
A check that survives only in the cold archive is flagged in the notes column.

Loads (front matter, worktree): literature-researcher `research`, omitClaudeMd; numerics-investigator instrument-doctrine,
nexus-debugging, nexus-impact, probe-cascade, vv-principles, numerical-bug-signatures, coding-elegance, rules loaded; explorer
nexus-exploring, nexus-guide, retirement-audit, omitClaudeMd; archivist retirement-audit, instrument-doctrine, nexus-verification,
nexus-exploring, vv-principles, algebra-of-record, doc-issue-relocation, rules loaded; method-implementer retirement-audit,
vv-principles, cross-domain-frames, algebra-of-record, coding-elegance, rules loaded; cross-domain-attacker cross-domain-frames,
nexus-exploring, omitClaudeMd; elegance-enforcer retirement-audit, coding-elegance, nexus-elegance, rules loaded; test-architect
retirement-audit, instrument-doctrine, nexus-verification, vv-principles, coding-elegance, algebra-of-record, rules loaded.

## literature-researcher — 21 of 21 COVERED rows read

Retired bodies: `git show 21f999c0:.claude/agent-memory/literature-researcher/lessons.md` (542 lines; at `1039683f` the file is
an 11-line retirement note). Carriers: `AGENT.md` (definition), `skills/research/SKILL.md`. Positive control: `Phantom`, 1 hit in
the skill.

| row | verdict | evidence / lost check |
|---|---|---|
| L-001 base | CARRIED | definition §3 "name the geometry in every mapping"; the routing advice (start at BMC for the ΔA/w factor, not Lewis-Miller) is in `bailey_ane_chimera_citation_refuted.md:97` |
| L-001 correction | CARRIED | definition §3 "Map a direction cosine by the operator it multiplies, never by its letter" is the redistribution-term test; the table is corrected in §3 |
| L-002 base | CARRIED | definition §2 Provenance; skill "Phantom-citation guard" |
| L-002 sharpening | CARRIED | definition §2 (slot, volume against year, "find the paper its equation numbers fit"); steps (a), (c) kept in new L-002 |
| L-003 base | **PARTIAL** | see below |
| L-004 | CARRIED | definition §2 "Method, not topic"; skill "Catalogue ≠ method-source" |
| L-005 base | CARRIED | skill "classify a method by READING its derivation, not by which catalogue cites it" (the definition alone would not carry it: its "Method, not topic" body speaks only of a catalogue and a geometry) |
| L-005 sharpening | CARRIED | definition §2 "the paper's own transport equation settles it"; the future-work clause kept in new L-013 |
| L-006 | CARRIED | definition §1 and "Not in the local library" ("never a pivot") |
| L-007 | CARRIED | definition §1 item 2 (refused connection, the memo note) |
| L-008 | CARRIED | definition §1 item 1 (`mocr.json` path; "Search the sidecar; verify on the rendered page") |
| L-009 base | CARRIED | definition §1 (`oa_status`); status values and the moved archive kept in new L-009 |
| L-010 base | CARRIED | definition §2 "The equation" (twin, seed or terminus, table, report the typo); vacuity kept in new L-010 |
| L-012 base (a) | CARRIED (by pointer) | the rule is not loaded; new L-012 preamble keeps "The count of sources behind a negative is the instrument-doctrine rule's statement X2". Note: a pointer sentence, not a `check:` line; see the optional restoration below |
| L-012 sharpening 08-12 | CARRIED | definition §1 "Searching" (probe a known paper, homographs, ISSN scope, citation graph) |
| L-012 second case (a) | CARRIED | definition Return "list those sites"; (b) corpus growth kept in new L-012 |
| L-013 base (a),(b) | CARRIED | definition §2 "A named scheme"; (c) and the rank check kept in new L-013 |
| L-014 | CARRIED | definition §2 "A disagreement between two sources" (edge, trade, missing edge = open) |
| L-015 (a) | CARRIED | definition opening paragraph and §1 |
| L-015 (b) | CARRIED | definition §1 item 1 (Python, NUL byte); page-marker recipe in `endf6_gendf_njoy_n2n_formats.md:44` |
| L-015 rider | CARRIED | definition §2 "even in a born-digital PDF, where layout loss…" |

**L-003 base, the lost check.** Retired text: "when a code comment names a reference, confirm the equation actually lives in THAT
paper before reusing the citation in output. Flag the mismatch to the user." The founding case cited a REAL paper (Adams-Yang-Zika
2008 JCP, a diffusion paper) for machinery that is in BMC 2010. The definition's Provenance clause is satisfied once "a database
confirms it", which a real-but-wrong paper passes; "find the paper its equation numbers fit" applies only to a citation already
known false; new L-003 checks a section pointer only "when a brief hands over a section and an equation number together".
Predicate re-run over the definition, the skill directory and the whole memory directory (pattern
`lives in (that|the cited)|wrong paper|equation actually|equation is in|in the cited paper|cited paper`; control `database
confirms`, 1 hit): 1 hit, `sood_fn_method_full_extraction.md:157`, another sense.
Restore to the head of new L-003:

> - check: a citation that resolves in a database can still be the wrong paper (a real diffusion paper cited for angular SN
>   machinery). Before reusing a citation taken from a docstring or a brief, open the cited paper and find the equation or the
>   construction in it; when it is not there, report the mismatch with the sites that carry it.

Optional (L-012 base (a), not a defect): turn the preamble pointer into a check line, "check: a negative states its count,
`3 of 3 local sources that define a face cosine`, never a bare `no source`."

## numerics-investigator — 15 of 15 COVERED rows read (stratum A); sub-points of TIGHTEN rows: see the sample section

Retired bodies: M1, M2, M5 from `1039683f`; M3, M4, M6, M7, M8, L9 from `54809381^` (the parent of the 2026-09-22 retirement,
710 lines). L1, L2, L6, L11, L12, L13 were already "DUPLICATE" stubs at both commits (the 2026-09-21 distillation's own carrier
claim), so their bodies were read one distillation earlier, at `02e0ab1c^` (1228 lines): every check there was located too
(L2's "confirm streaming alone at α = 0, then the per-ordinate flat-flux test" is definition §3 steps 4 and 6 with Signature 1's
probe; L6's "Krylov-on-apply converges to the buggy operator's self-consistent fixed point" is the reason for step 6; L11's
c-dependent-error trigger is Signature 9; L12's three probes are `vv-principles` § Bit-identity and #5; L13's copy-not-checkout
isolation is Signature 11's probe). Kept-entry claims re-checked in the new digest: L14 (line 164, seed-blind ladder line 178),
L15 with `restart` (181), L16 (193), L17 (211), L21 (305), L23 (363), L25 (450). Carriers: definition (§1 to §4, Return), the
seven preloaded skills, the always-on rules. New memory: the worktree `lessons.md` "Retired IDs" stubs and kept entries.
Positive control: `OWN control`, 1 hit in new `lessons.md:532` (L26·8).

| row | verdict | evidence |
|---|---|---|
| M1 | CARRIED | definition §2 (rate, contract, kernel, ownership, ordering, transfer gain) and §4 "build the difference operator" |
| M2 | CARRIED | definition §3 step 7 (N at fixed fine mesh; h → 0); `vv-principles` Mode 7; the seed-blind MMS ladder instance kept in L14 |
| M5 | CARRIED | definition §3 "every probe budgets for its positive control"; instrument-doctrine X1; "each structural zero needs its OWN control" kept at L26·8 (`lessons.md:532`) |
| M3 | CARRIED | definition §3 step 1 (record against budget before a tolerance sweep), §2 hang; L10c (budget is dimension-dependent) and L15(b) (`restart`) kept |
| M4 | CARRIED | `vv-principles` Mode 7 (activates / nulls), Mode 12 (space-side stabiliser by degenerating the fixture); instances kept in L7·1, L8, L17·2, L23·8 |
| M6 | CARRIED | definition §4 "Two of your measurements contradict … re-measure on the subspace the driver actually feeds" |
| M7 | CARRIED | definition §3 step 1 (paper's approximation level); instances L4, L21·5–6, L25·5–6 kept; the docstring-as-claim instance is instrument-doctrine X3 |
| M8 | CARRIED | definition §3 step 1 (contaminated reference), §3 opening (mixture consistency identity, L20·5's content); `vv-principles` #6 (doc plus docstring are one source); L3, L16·3, L22 kept |
| L1 | CARRIED | `probe-cascade` "Don't skip probes" (`SKILL.md:126`); the "six wrong hypotheses" figure is war story, not a check |
| L2 | CARRIED | `numerical-bug-signatures` Signature 1 |
| L6 | CARRIED | definition §3 step 6 ("against a hand reference on a non-flat ψ"); record file named in the stub |
| L9 | CARRIED | definition §2 "A hang or a timeout"; the uncarried clause (reproduce on a fast proxy and bisect) is restored in the stub |
| L11 | CARRIED | definition §3 step 3; Signature 9; the API paths kept in the stub |
| L12 | CARRIED | `vv-principles` #5 second check (order at the cell; a cell at target order with the term off cannot be repaired) and § Bit-identity (swap, term-silent, AMPLIFY) |
| L13 | CARRIED | `numerical-bug-signatures` Signature 11 |

## explorer — 12 of 12 COVERED rows read (11 lessons rows plus the `MEMORY.md` section-2 row; the owner's "11" counts lessons only)

Retired bodies: L-001, L-002, L-004, L-005, L-007, L-039, L-046 from `lessons_archive.md` at `1039683f`; M-1, M-4, M-5, M-6 from
`ad16622a^:.claude/agent-memory/explorer/lessons.md` (the 2026-09-22 parent; at `1039683f` only M-3 survives). Carriers: definition
(§1 to §5, Return) and the three preloaded skills (`nexus-exploring`, `nexus-guide`, `retirement-audit`); no rule body (omitClaudeMd).
Positive control: `conflict to adjudicate`, 1 hit in new `lessons.md:267`.

| row | verdict | evidence |
|---|---|---|
| L-001 | CARRIED | definition §4 four searches; `retirement-audit` A.1 tells (`cached_property`-reached leaf, class-name bypass consumer), A.3 constructors; §2 Protocol-typed receivers and grep as primary evidence; the `singledispatch` arm half has no object: `[M]` AST pass over 404 files of `orpheus/`, 0 decorators containing `singledispatch` or `.register` (control: 1322 `property`/`dataclass` decorators) |
| L-002 | CARRIED | definition §3 "an issue's premise before planning its work … close-and-verify" |
| L-004 | CARRIED | definition Return "A verdict names both arms, the discriminator … the value judgement goes to the orchestrator" |
| L-005 | CARRIED | definition §3 `merge-base --is-ancestor` |
| L-007 | CARRIED | (1) re-run the census last and (3) mark the line numbers' vintage: new L-046 (`lessons.md:169-176`, "re-derive every cited `file:line` from its anchor last … mark each successor HEAD vs working tree") with Return's anchor-string rule; (2) diff the uncommitted hunks: new L-012 (`lessons.md:270`) |
| L-039 | CARRIED | definition §2 "its `unresolved` count is the census", "grep or an AST pass is primary evidence"; two-level resolver is `retirement-audit` HELPER-RETURNED-RECEIVER |
| L-046 (archive, degree) | CARRIED | definition §4 "A hub is not a template … read the convention it encodes against the derivation of record" |
| M-1 | CARRIED | definition §3 "cheapest decisive probe"; instances kept (L-011 "handles it", L-020 negatives, L-036 plan staleness); the test self-description instance is `retirement-audit` D.14 ("are the two sides still INDEPENDENTLY produced?") |
| M-4 | CARRIED | definition §2 last paragraph (fields, methods, callables, Protocol receivers, captured functions, properties, docstring roles); the caller-module binding kept in new `lessons.md:167` |
| M-5 | CARRIED | definition §2 last bullet (run, control, inert run); the `None` arm and same-data rebuild kept (`lessons.md:223`, L-041) |
| M-6 | CARRIED | definition Return; "conflict to adjudicate" kept in L-036 (`lessons.md:267`) |
| `MEMORY.md` §2 active state | CARRIED | definition §3 `is-ancestor`; the remainder is campaign narrative |

## archivist — 8 of 8 COVERED rows read

Retired bodies: `1039683f:.claude/agent-memory/archivist/lessons.md` (572 lines; SPINE line 13, §1b, §2b, §6 line 463, §7, §9).
Carriers: definition, the seven preloaded skills, the always-on rules (path-scoped `documentation` loads for `docs/theory/**`).
Positive control: `forbidden sentence`, kept at new `lessons.md:435`.

| row | verdict | evidence |
|---|---|---|
| THE SPINE | CARRIED | definition §1 (live tree this session; V&V vocabulary verbatim), §2 (warning set, never a count; dead references), §4 retirement by tense |
| 1b.5 negatives re-verified at END | CARRIED | definition Return "every negative you are about to publish ('not landed', 'no consumer', 'no pin')" |
| 2b.1 xref gate role-blind | CARRIED | definition §2 "certifies `:mod:` targets only until #497 lands: accept a page on your own import probe with a live and a dead control"; mechanism text in archive L-125 |
| 6.23 architecture-pass ladder | CARRIED | definition §1 "the production docstrings, the diagnostic scripts (read and run them), the commit bodies and the design memos, in that order. Numerical bounds come from the test files" |
| §7 preamble | CARRIED | definition §1 V&V vocabulary verbatim; the forbidden-sentence bullet stays in the digest |
| 9.2 `_build` filter | CARRIED | `code-search` (`_build` exclusion), definition §2 `git grep`; the "check Y is there" tail is definition §1 "checked against the live tree … not the brief" |
| 9.3 error-entry without catches | CARRIED | definition §5, near-verbatim |
| 9.5 raise string executable | CARRIED | definition §5 ("a string inside a `raise` … is code: report it"); `retirement-audit` B.11 shortest distinctive fragment |

## method-implementer — 6 of 6 COVERED rows read

Retired bodies: `2990cebb^:.claude/agent-memory/method-implementer/lessons.md` (424 lines). Carriers: definition (Done), the five
preloaded skills, the always-on rules. Positive control: `weights @ stacked`, 1 hit in `coding-elegance`.

| row | verdict | evidence / lost check |
|---|---|---|
| L-002 pyright scope traps | CARRIED | definition Done "`npx pyright <touched files>` … against its `total: 0` baseline" measures the seam file, which also catches debt imported by a promotion into a scanned directory |
| L-003 orphaned sole-use imports | **PARTIAL** | see below |
| L-004 three searches | CARRIED | `retirement-audit` A.1 tell (class-name bypass consumers, direct constructors), A.3 |
| L-005 never `git checkout` a mutated WIP file | CARRIED | `process-discipline` "Mutation-testing an uncommitted file" |
| L-007 sub-floor term, structural teeth | CARRIED | `vv-principles` Mode 10 (machine-precision threading, consumed sign flip, no-op control, "NO value-improvement leg — do not manufacture one") |
| L-012 bare-ndarray shape lies | CARRIED | `coding-elegance` anti-pattern 20 special case "verify the ACTUAL shape at the seam" |

**L-003, the lost check.** Retired text: "after deleting a symbol, grep the module for the short-name of every import the deleted
body used (`spla.`, `sp.`, the helper name). A pyright `reportUnusedImport` will also catch it". The owner's carrier is
`retirement-audit` G.24 read in paraphrase; G.24's instrument is a dependency-audit table of helper symbols written before the
plan, and its tell is "an orphaned helper", so an import binding orphaned by the delete is outside what it asks for. The
retired text's backstop is doubtful in this tree: `[tool.pyright]` in `pyproject.toml` sets neither `reportUnusedImport` nor
`typeCheckingMode` (`[M]` read 2026-11-16), and pyright's default leaves that rule off (`[R]`, not probed). Predicate re-run over
the definition, the five preloaded skills, every rule (`coding-standards` included) and the whole memory directory (pattern
`unused import|reportUnusedImport|sole-use|orphan(ed)? import|import.{0,30}(orphan|unused)`; control `orphaned helper`, 1 hit in
`retirement-audit`): 4 hits, all in two topic files recording the founding instance (`issue_257_s7_scipy_single_source_closeout.md:134`,
`:150`; `issue_240_d5b_s3_unified_matvec_closeout.md:231`), none stating the check. Severity: low (dead code, no wrong answer).
Restore to the digest:

> - check: after deleting a symbol, grep its module for the short name of every import the deleted body used (`spla.`, `sp.`);
>   an import with no remaining use is part of the retirement (`retirement-audit` G.24).

## cross-domain-attacker — 5 of 5 COVERED rows read (4 plus item 10 "in part")

Retired bodies: `69bfc44e^:.claude/agent-memory/cross-domain-attacker/lessons.md` (411 lines; M1 to M4 at lines 25-194).
Carriers: definition (steps 0 to 5, two standing bars, Output), `cross-domain-frames` `SKILL.md`, and `reference.md`, which is not
preloaded but which the skill and definition steps 1 and 2 direct the agent to (Part A, Part C); no rule body (omitClaudeMd).
Positive control: `two-question`, 2 hits in `reference.md`.

| row | verdict | evidence / archive-only riders |
|---|---|---|
| M1 items 1, 2, 3, 5 | CARRIED, riders archive-only | definition step 0 (branches, fibre, consumers, fused attributes, straddler inversions, optimality theorem). Four riders live only in `lessons_archive.md`, to which the digest points by number (L-017, L-022, L-023): an existing satellite `cached_property` is the strongest counter-evidence to a co-sourcing factory (archive 1120); grep `<b>\.<a>\._` and `setattr` before accepting a save-state charter (1217); report the row you got wrong by following a chartered rule literally (1203); chase the optimality theorem's equality (alignment) condition (815) |
| M1 item 6 | CARRIED | `reference.md` Part C Smell #17, shapes (a) and (b), read at the clause |
| M1 item 10 | CARRIED | step 5 (word already spent); the verb grep restored into the digest's M1 paragraph; the diagonal-multiplier companion restored in the Part 3 multiplier row (`lessons.md:500`) |
| M2, M3, M4 | CARRIED, riders archive-only | first standing bar (both directions, RHS spellable, the violating object), second bar (refuted FOR a question), Output ("a call site, a widened return type, a verb on an existing type, or a gate"). Archive-only riders: ask what the chartered hypothesis is necessary and sufficient FOR, since naming the adjacent law is the deliverable (854); a guarded normalisation constant usually CANCELS in the law (877); which spelling of each face the construction uses (873); the steelman at the smaller scale by currying (1142); zero production callers moves the question to the consumer side (1129); a stale-owner docstring is a free finding (1144) |
| P4 affine/torsors | CARRIED (rests on `reference.md`) | `reference.md` A.1 torsor row (overturned example, two-question test); the new digest keeps no torsor line (`two-question` hits only the archive and `reference.md`), so this verdict holds only if `reference.md` counts as a carrier (not preloaded; the skill's `SKILL.md` lines 22, 37, 62 send the agent to it) |

No row is PARTIAL by the brief's definition, since the archive is a file of the owner's memory and the digest names its
L-numbers. The flag stands because the archive is cold: of the ten archive-only riders, the three that change a verdict rather
than add colour are the necessary-and-sufficient question (M2), report-the-row-you-got-wrong (M1.3) and zero-production-callers
(M4). Optional restoration, one line each, into the digest's M1 paragraph:

> Riders that change a verdict (archive): ask what a chartered hypothesis is necessary and sufficient FOR, and name that adjacent
> law (→ L-018); report the row a chartered rule made you get wrong (→ L-023(b)); shipped machinery with zero production callers
> moves the question to why nothing calls it (→ L-022).

## elegance-enforcer — 3 of 3 COVERED rows read

Retired bodies: `1039683f:.claude/agent-memory/elegance-enforcer/lessons.md`, "Standing review order" 0 to 4 (lines 22-38).
Carriers: definition, preloaded `retirement-audit`, `coding-elegance`, `nexus-elegance`, the always-on rules. Positive control:
`reassur`, 1 hit in `rules/plan-authoring.md` (another sense) and 1 in `nexus_architecture_survey.md:59`.

| row | verdict | evidence / lost check |
|---|---|---|
| Standing order 0 (adversarial first) | **PARTIAL** | see below |
| Standing order 3 (filter `_build`) | CARRIED | `code-search` (`docs/_build` exclusion), definition §1 "excluding `docs/_build`" |
| Standing order 4 (three legs) | CARRIED | definition §4 "A VIOLATION needs all three; with one missing it is a CONCERN" |

**Standing order 0, the lost check.** Retired text: "Ask 'how would I BREAK this' (silent wrongness, weighted to defects that
fail in the REASSURING direction) and 'how would I make this 100× better' (reframe the JOB)." Definition §3 and the workflows
rule carry the two questions and the withdrawn-attack form, but not the weighting: "the reassuring direction first" is in
`docs/development/workflows.md:130` (linked from the definition, not loaded) and in qa's definition, and in this owner's memory only
in `nexus_architecture_survey.md:59`, whose new index hook (7 words) says the file's method lesson is "def §3 (COVERED), not a
reason to open the file", so nothing routes to it.

The retired line pointed at "L-020", which is absent from the old digest, the archive and the new memory (`grep -rn L-020` over
the owner's directory: 0). `git log -S 'L-020'` finds it retired on 2026-09-22 by `8d80dc29`; its body
(`git show 8d80dc29^:.claude/agent-memory/elegance-enforcer/lessons.md`, "### L-020 — A balanced survey is the SECOND phase")
carries five checks, of which definition §3 and the workflows rule (W5, "a second-pass verdict never bounds the first") carry the
two-pass order, the reframing and the no-bounding rule. Not carried anywhere the owner loads or in its new memory (`[M]` grep of
the definition and memory for `object-shaped|what objects are hiding|ONE absence|no live instance|result type`: 0 relevant):
(a) weight the reassuring direction ("a short list, an empty result, a `None` — because the caller cannot tell it from
correct"); (b) "when a codebase STORES `confidence`/`source` on its data, grep whether the RESULT type carries them"; (c) the 100×
answer "explains several independent break-it findings as ONE absence"; (d) separate "this was by design and I hadn't read it" (a
withdrawal) from "the mechanism is real but no live instance exists" (a surviving CONCERN); (e) ask "what objects are hiding here,
and what context do they share?", not "which files do these methods go into?". So this row is the 2026-09-22 false-retirement
class, one step removed. Restore L-020 to the digest from `8d80dc29^` in its own words (war story to the archive), at minimum:

> - **L-020** (restored 2026-11-16). The break-it question is "what input makes this confidently return a WRONG answer?";
>   weight the reassuring direction (a short list, an empty result, a `None`), because the caller cannot tell it from correct.
>   When the data stores `confidence` or `source`, grep whether the RESULT type carries them. A 100× answer explains several
>   independent break-it findings as ONE absence. A withdrawal ("by design, I had not read it") differs from a surviving CONCERN
>   ("the mechanism is real, no live instance exists"). Ask what objects are hiding and what context they share, not which files
>   the methods go into.

and UPLIFT (a) into the definition §3 source, `docs/development/agents/elegance-enforcer.md`, in qa's words: "(the input that
makes it silently wrong, the reassuring direction first)".

## test-architect — 1 of 1 COVERED row read

Retired body: `1039683f:.claude/agent-memory/test-architect/lessons.md:18-21`.

| row | verdict | evidence |
|---|---|---|
| H.2 THE SPINE | CARRIED | definition §5 "Every gate names the input in today's tree that reddens it … under `python -O -m pytest`" (a), §3 "REFERENCE (a structurally independent route)" (b), §3 "Activation" (c) |

## Stratum B — COVERED sub-points inside TIGHTEN rows: sampled, 11 of about 35

Not a census: about 35 sub-points are marked COVERED inside TIGHTEN rows (numerics-investigator about 27 across L7, L10, L15,
L16, L18–L28; test-architect 3; cross-domain-attacker 1; method-implementer 1; a few elsewhere). The sample is biased toward
rows whose carrier sentence the owner quoted; the unsampled "core COVERED" shape (a multi-point lesson mapped to one definition
sentence: NI L18 four clauses, L19·1, L23 ·2/·9/·10, L24·1, L25·1, L26·1) is the stratum most like the 2026-09-22 failures and
is unread. Sampled, largest owner first:

| owner, sub-point | verdict | carrier |
|---|---|---|
| NI L22·1 | CARRIED | `vv-principles` § Bit-identity "print `max|a−b| / max|b|` BEFORE reading the nulp count" |
| NI L22·3 | CARRIED | same section, "a producer's TOLERANCE pin can never warn its consumers' BIT-IDENTITY pins … byte-level fingerprint" |
| NI L22·5 | CARRIED | `code-search` (`--color=no`; `sys.path[0]`; print `module.__file__` first) |
| NI L28·2 | CARRIED | `vv-principles` #38 "±3-ULP jitter over ≥10 seeds; a flipped verdict is noise" |
| NI L7·2, L7·3 | CARRIED | definition §3 step 4 (independent kernel bit for bit; convention at a producer–consumer seam); stubs kept |
| TA 2.11 traced-memo half | CARRIED | definition §5 "bypass() at pytest_configure, or a textual arm"; the spy-the-consumer's-binding half kept |
| TA 5.1 absolute-atol half | CARRIED | definition §3 "a residual is normalised by what it divides by" |
| TA 7.5 `replace` half | CARRIED | `coding-elegance` Pattern 4 corollary (`dataclasses.replace` re-runs `__post_init__`) |
| CDA D2(c) | CARRIED | definition step 4 "A frame keyed to an operator's algebraic shape fires only on the members with that shape" |
| MI L-001 bullets 1–3 | CARRIED | `coding-elegance` anti-pattern 19 (variance, `TYPE_CHECKING`, `weights @ stacked`) |

## Totals

| owner | read | CARRIED | PARTIAL | FALSE | archive-only flag |
|---|---|---|---|---|---|
| literature-researcher | 21 of 21 | 20 | 1 (L-003 base) | 0 | 0 |
| numerics-investigator | 15 of 15 | 15 | 0 | 0 | 0 |
| explorer | 12 of 12 | 12 | 0 | 0 | 0 |
| archivist | 8 of 8 | 8 | 0 | 0 | 0 |
| method-implementer | 6 of 6 | 5 | 1 (L-003) | 0 | 0 |
| cross-domain-attacker | 5 of 5 | 5 | 0 | 0 | 2 rows (10 riders) |
| elegance-enforcer | 3 of 3 | 2 | 1 (standing order 0) | 0 | 0 |
| test-architect | 1 of 1 | 1 | 0 | 0 | 0 |
| stratum B sample | 11 of ~30 | 11 | 0 | 0 | 0 |

71 of 71 fully COVERED rows read at their bodies: 68 CARRIED, 3 PARTIAL, 0 FALSE. Unlike the 2026-09-22 retirements (false in part for five
owners), today's COVERED claims hold at the clause; the three misses are each one check with a carrier that names the topic
but not the instrument.

Stamps: worktree HEAD `bcd4478b` at start and at end; worktree clean at both (`git status --short` empty). Nothing tracked was
edited; this file is the only write.
