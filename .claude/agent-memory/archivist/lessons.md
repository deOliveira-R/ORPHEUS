# Archivist — Lessons (hot digest)

Read FIRST, every dispatch. One imperative per lesson with the failure→correction core that earned
it, and a `→ L-0NN` pointer into **`lessons_archive.md`**, which holds the war story, the `[M]`
numbers and the `file:line` detail — open a section on demand, never the file.

**Three things are NOT repeated here; a line that restates one is retired on sight.** (1) A rule's
or skill's clause — `instrument-doctrine` X1–X4, `articulation`, `process-discipline`,
`code-search`, `retirement-audit`, `vv-principles`, `plan-authoring`, `coding-standards`,
`vv-testing` — cited by ID, never quoted. (2) Mechanical HOW: build-gating, venv/worktree facts and
the 9-step close-out arc are `AGENT.md`. (3) Branch-1/Branch-2 discipline: `algebra-of-record`.

**THE SPINE — a page is DONE when** every cross-ref resolves against the LIVE tree · every claim
was verified against live code THIS session · every claim's V&V level matches `vv-principles`
verbatim · every retired symbol leaves no present-tense-false mention · the build's
WARNING/ERROR/CRITICAL **set** is unchanged from a freshly-measured `-E` baseline. Every rule below
is one face of that; AGENT.md's Quality Checklist is the same bar, itemised.

---

## 1. Ground truth is the LIVE tree — every other surface lies eventually

**The meta-rule (AGENT.md checklist #6): the brief is the FLOOR, live code is the rule.** Brief,
docstring, verdict memo, design review, commit body, retirement shim, scanner finding, plan row,
"MEASURED — do not re-derive" block and a mid-task delta are point-in-time snapshots. Verify, then
write, then FLAG every scope-expansion the verification forced. All of §1 is a way that fails.

### 1a. Reading a relayed claim

- **A brief can instruct a HARD audit error — read the scanner before obeying a `vv-status`
  instruction.** `documented` is the ONLY legal status (`tests/_harness/audit.py`, exit 2
  otherwise); what a new `verifies` marker earns is UN-SENTINELING, not an upgrade. → L-106
- **A brief describing code being typed NOW is a FORECAST** — verify the symbol, and poll on the
  INVARIANT (`! git grep -q <old> -- orpheus/`), never on one new symbol: a half-written tree
  answers plausibly, and its `AttributeError` reads as "unsupported", not "broken". → L-087, L-098
### 1b. The tree moves under you

- **On a live branch the re-read is a LOOP, not a pre-flight — re-read the public surface after
  EVERY build, and run `git log -3` + `git status --porcelain -- orpheus/` at the END.** A clean
  `orpheus/` late means the concurrent pass COMMITTED and its message is the diff of your premises.
  Highest-decay class: a field's `compare`/`repr`/default, a guard's clause COUNT, a helper's
  SIGNATURE — `-W` is silent on all three. → L-089, L-088, L-094, L-097, L-085
- **A module under concurrent edit: publish NOTHING you count in it** — point at the artefact that
  re-measures it; if a count must be stated, date it and name the census. → L-095, L-088
- **Discriminate your own edits from a concurrent editor's by CONTENT, never by the porcelain flag
  or a shared date.** → L-070, L-056
- **A gap YOU report has the shortest shelf life in the corpus — its own commit can close it**,
  because comparing two implementations is simultaneously what exposes a gap and what motivates the
  repair. Re-run every gap claim after the session's LAST code edit and publish it as
  history-with-its-repair-hash, which cannot rot. → L-080, L-075, L-087
- **Every "NOT landed / still open / zero consumers / no pin / untestable / unreachable" clause is a
  negative claim about a tree someone else is editing — re-verify at the END.** → L-107, L-073,
  L-082, L-079

### 1c. Numbers you publish

### 1d. Read the object, not its description

- **Construct the object and print the property.** `hasattr` is False for a dataclass field with no
  class default and for an `__init__`-assigned attribute; the ladder is `hasattr` →
  `dataclasses.fields`/`__annotations__` across `__mro__` → `self.x =` across the MRO → CONSTRUCT.
  → L-093, L-076, L-053, L-105
- **For "all N of family F do X", re-derive F's MEMBERS, not |F|** — a family turns over with its
  count fixed, and a roster's own list is not the family (`vars(cls)` +
  `isinstance(v, classmethod)` said FIVE where four pages said four; the fifth row was the
  argument). → L-101, L-074, L-081
- **The pre-carve tree is measurable** (`git worktree add … HEAD --detach`, or `git archive HEAD`),
  and when the retired object is shape-characterisable a hand-built STAND-IN beats a diff: the
  equal-shape CONTROL is the finding, and a diff cannot produce it. Print `orpheus.__file__` as proof:
  `PYTHONPATH` beats the editable install, but a probe run from inside the repository loses to
  `sys.path[0]` (`[REFUTED 2026-09-22]` the earlier claim that the editable install outranks
  `PYTHONPATH`; the mechanism is `code-search`'s). → L-095, L-101, L-050, L-083
- **Any cross-carve table owes a line naming WHICH statistic is comparable** — the tell that you
  need it is that the CONTROL column moved. → L-095

### 1e. Surfaces nothing gates

- **A published `.. code-block:: python` is the highest-severity staleness there is** — it promises
  reproducibility and nothing gates it. After any constructor/signature change grep the code-block
  BODIES before the prose sweep, and RUN what a present-tense sentence promises works. → L-077,
  L-047, L-093, L-099, L-105
- **A regenerated cache or untracked artefact is a doc surface with no symbol in it; `ls -l` is the
  only instrument** — `git status`, `-W`, the xref gate and `dead_references` are all blind. → L-093
- **A page's own ⚠ CAVEAT can BE the diagnosis a later carve lands, or the theorem the next step is
  built on — keep it verbatim under a ⭐, never tense-flip it.** Tell: *"X is not always Y, so the
  precise statement is Z, not a named attribute"* is a page saying where its subject is welded.
  → L-104, L-086

### 1f. Landings, predictions and seams

- **When a fix OVERRIDES an artefact that still exists and will be read FIRST, the page is a
  COUNTER-RECORD, not a report.** A literature memo's implementation note or a plan's pseudocode
  survives the fix and outranks it in a future search, so give the correction its own named
  subsection and carry the falsifying table VERBATIM (the no-op row, the pass, the redundant, the
  degenerate-coincidence). Without the table the next session re-applies the wrong injection
  point, because the surviving memo says so. The doc-side of `process-discipline`'s
  "a refuted candidate is first-class output": the structural reason must be findable *before*
  the artefact it refutes. → L-114
- **A two-clause ruling may have an occupant for only one clause** — publish clause 1 as the ruling
  and clause 2 with the census showing it empty. → L-108
### 1g. Instruments and corroboration

## 2. The build is BLIND to most doc-correctness defects — grep is the gate

**Meta-rule: `-W` proves only "I added no NEW warning". Acceptance for a correctness sweep is a
grep inventory with a per-hit KEEP/FIX adjudication.**

### 2a. What is silent, and what actually warns

- **Silent at DEFAULT severity:** an unresolvable `:func:`/`:class:`/`:meth:`/`:attr:`/`:mod:`; a
  role into a `:noindex:` or un-`automodule`'d module; every violation of RST's
  no-nested-inline-markup rule (`**``lit``**` renders its delimiters; a role opening straight after
  `. * ~ § ↔ =` dies outright, eating the LaTeX backslash); a trailing space before a closing role
  backtick (it swallows the sentence); a backslash inside a `` ``literal`` `` (scientific notation
  is a `:math:` role, never a literal); a raw relative hyperlink; a raw file path; a kwarg or
  expression inside a `code-block`. → L-002, L-044, L-069, L-074, L-082, L-058
- **These DO warn, so they are gated:** a dangling `:eq:` (measured with controls) and — per the
  refutation above — a dangling `:ref:` and `:doc:`, cross-doc included, so BOTH label classes are
  rename-gated by the build; a bad `:by:`; an unresolved
  `catches("ERR-NNN")`; `ref.ref` *"A title or caption not found"* (a bare `:ref:` to an anchor
  above a paragraph, an admonition, or a **bold run-in heading**); a `:widths:`/column mismatch; a
  malformed `===` table (never hand-align one holding a `:math:` role — source length is not
  rendered length); an italic run interrupted by a role; **an ORPHANED citation** — a `[Key]_`
  defined but no longer referenced warns `Citation [Key] is not referenced. [ref.citation]`, so
  the DELETION direction is gated where §9's `grep '^\.\. \[Key\]'` covers only the adding one
  (a cut that removes the last citing section moves the citation to a surviving section on the
  same topic, or deletes it when a sibling already carries its own); **and a TRAILING UNDERSCORE
  in prose** — `Γ_-`, `S_-`, `X_+` are RST reference syntax and raise `ERROR: Unknown target
  name`, on which `-W` exits 1 (the escaped `Γ\_-` and the no-trailing-`_` `V_bulk` are silent;
  `[M]` all 7 corpus occurrences already sit inside a literal, a comment or a directive option).
  ⚠ `SyntaxWarning` and the nexus per-marker
  lines carry no `WARNING:` prefix, so read the whole log and grep `CRITICAL:` too. → L-114,
  L-070, L-060, L-002, L-027, L-040, L-054, L-055, L-048, L-095

### 2b. The xref gate, and what acceptance actually is

- **⛔⛔ The project xref gate is ROLE-scoped-blind: `DEAD TARGETS: 0` certifies `:mod:` targets and
  NOTHING else, in `.rst` and in `.py` docstrings alike.** `judge()` re-checks the target's HEAD
  carrying the ORIGINAL role, so every dead fully-qualified `orpheus.*` target under a non-`mod`
  role reads DECLINED. `[M]` 2026-09-17: a dead `:class:`/`:attr:`/`:meth:` → DECLINED while a dead
  `mod` → DEAD, and the gate printed `DEAD TARGETS: 0` over four live dead `:class:` roles. **The
  2026-08-24 "appears repaired" note is REFUTED** — it rested on two instruments sharing the
  blindness. ⟹ acceptance for a page is YOUR OWN import probe carrying BOTH controls (a live target
  must read ALIVE, the retired one DEAD); a bare `0 dead` is unreadable without the negative
  control. The one-line fix (`head_role = "mod" if "." in target else role`) is still UNLANDED.
  → L-111, L-082, L-067, L-062
- **What the gate does cover, and its false negatives:** it resolves FULLY-QUALIFIED roles by
  IMPORT (render coverage irrelevant), reads whole `.rst`, walks only `orpheus tests docs`, judges
  only `DECIDABLE_ROOTS`, skips UNQUALIFIED roles by design, and is blind to LITERALS. A PEP-420
  namespace package imports fine with 0 members; a role wrapped INSIDE its dotted path is skipped
  AND renders plain text. Never touch its ALLOWLIST. Measure a patch as a COPY at depth 1 inside
  the repo, run as a SUBPROCESS, with a throwaway `docs/_ctl.rst` control — **stock == patched is
  itself the tell that the patch is inert.** → L-045, L-046, L-047, L-052, L-062, L-071, L-082
- **Nexus `dead_references` resolves by RENDERED TARGET where the gate resolves by IMPORT; the SET
  DIFFERENCE is the triage** (hook-only ⇒ un-surfaced-but-live; both ⇒ really retired). It RESCUES
  an inherited member, which decides a sweep's width: a RENAMED/RETIRED member must re-point, a
  member that merely MOVED re-points only where the SENTENCE claims it is defined there. Before
  believing a dead target's NAME, read its graph EDGES (`type_uses` = a code bug, `documents` = a
  page, `references` = a docstring) — a name can be an artefact minted by a third tree. → L-094,
  L-052
### 2c. The render is the oracle for markup; the source diff is the free pre-check

- **The HTML slice is the only instrument that sees nested inline markup — and the shape YOU will
  write is `**``value``**` in a numeric `list-table` cell.** Strip tags → unescape → slice by
  `id=` with `rfind` (or `role="main"` … `<footer|sphinxsidebar`) → require 0 visible backticks and
  0 surviving `:role:` spellings. **Assert the slice contains known page prose AND check its
  LENGTH**, or a dead anchor reports "clean". Exclude `_modules/` and orphaned `_build` pages.
  → L-068, L-069, L-074, L-080, L-082, L-084, L-085
- **A SOURCE regex does gate it if it is `re.S` over CONTIGUOUS ADDED blocks and set-differenced
  against `git show HEAD:<file>`** — strip `code-block` bodies and literals first, cap the run
  length, and add a `**`/`*` PARITY check per new paragraph (skipping `.. math::`). ⛔ Do NOT rebuild
  the naive `\*\*(.+?)\*\*` "role inside bold" regex: it pairs one run's closing with the next's
  opening and floods (119–132 hits on clean files). → L-076, L-079, L-095, L-097, L-101, L-074
- **Keep the pre-edit `-E` build and DIFF it per page — a rendered delta needs no provenance
  argument.** On a tombstone-heavy page take the multiset DIFFERENCE and read `added == removed` as
  a CONTEXT SHIFT, not two events; a page-wide count indicts someone else's prose. ⚠ A WRAPPER's
  exit is not sphinx's, twice over: `nohup … &` inside a background Bash call reports the SHELL's,
  and `python -m sphinx -W … | tee log` reports **tee's**, so a failing build reads as a pass —
  redirect (`2>/tmp/log; echo $?`) rather than pipe. → L-072, L-088, L-089, L-114
- **`<cite>` in the built HTML is the Markdown-port smoking gun — but count both spellings before
  "fixing" it**: `` `[M]` `` is this corpus's marker and renders `<cite>`, so normalising yours
  makes your text the inconsistent one. A port's warning count is a non-representative sample of
  its defect count: census the delimiter alphabet (a 3+-backtick run-length histogram is a TOTAL
  census) and use the port's SOURCE as the oracle. → L-061, L-071
- **Smartquotes mis-directs a closing `"` after an inline literal — fix by EXTENDING the quoted
  fragment to end on a WORD.** And ⛔ tombstone a quoted claim with plain quotes + `, verbatim,`,
  never an outer `*…*`: bold or a literal nested inside an italic quotation leaks, `-W` silent.
  → L-082, L-085, L-086, L-056
- **For an un-`automodule`'d module the build sees nothing — substitute a DIFFERENTIAL docutils
  parse, HEAD vs working tree, counting roles that survive as TEXT.** Subtract the Sphinx-only
  classes first (`:label:` on `.. math::`, unknown roles) or the delta is all noise; bare docutils
  knows no Sphinx domain roles, so only `:math:` is testable this way. → L-078, L-105, L-096

---

## 3. A `:label:` is a V&V edge — grep the matrix before touching it

- **NEVER rename or delete a label a `@pytest.mark.verifies(...)` targets** — grep `orpheus/` and
  `tests/` FIRST; a hit means you report the edge (you do not edit `tests/`). The marker also
  decides WHICH BODY a label keeps: read the test body, let the existing label keep what its marker
  asserts, and mint the NEW label for the generalisation. → L-003, L-032, L-094
- **Four fates under an ontology overturn, decided by the label's BODY:** states the retired claim
  ⟹ retire + repoint every `:eq:` citer · still true ⟹ rename · untouched by the overturn ⟹ KEEP +
  a `.. note::` that the prefix is a historical artefact · a retired STEP whose CONCEPT survives ⟹
  keep-and-repoint. A stale NAME is not a false CLAIM; show a retired equation in the history
  section as an UNLABELLED `.. math::`. ⚠ The KEEP fate was argued partly from *"a cross-doc
  `:ref:` miss is silent"* — **REFUTED** (§2a): a dangling `:ref:` warns, so the rename risk is a
  build failure, not a silent break. KEEP still wins when the BODY is untouched (a stale name is
  not a false claim) and when the citers are many, but say so on those grounds. → L-063, L-003,
  L-076
- **A labelled equation stating a CALL SIGNATURE is a time bomb** — move the signature to prose and
  let the markers decide what the label keeps. Name the OBJECT, never the paper: an eq-label naming
  a citation is a latent staleness bug by construction. → L-108, L-070
- **Two labels for one equation: publish the REGISTER each page owns, do not collapse** —
  collapsing moves a generated matrix row and re-points markers a docs pass may not touch. → L-070,
  L-064
- **Classify every label you add**, and know the sentinel arithmetic: `documented` is a KIND, not a
  coverage statement (a label can honestly be verified by 9 tests AND sentineled), it must name a
  `:label:` in the SAME file, it works INDENTED, and it moves `matrix.rst`'s sentinel count without
  moving the test count. Algebra-of-record SymPy-identity labels are verifies-COVERED, not
  documented. Orphan adjudication: WIRE / SENTINEL (three shapes) / GAP — never manufacture a gap.
  → L-004, L-035, L-036, L-039, L-065, L-077, L-063, L-049, L-037, L-027, L-030
- **Declaring `.. implements::` switches token-inference OFF for the WHOLE equation, so an
  incomplete declaration UNDER-covers** — ask *what else computes this?* before the first directive
  (7 of 14 needed 2–4). Writing the explanation MINTS new guesses, so never publish a live guess
  count. A dead `:by:` has four fates (migrate · remove · remove-and-add-the-new-home · keep BOTH
  as declared delegations); a RED `-E` baseline's `nexus.directive` warnings ARE the retirement
  site list; predict and check the `wrote N edges` delta. An equation with no implementer keeps its
  guesses forever, so `no-implementation` is a KIND worth a section (identity · law ·
  canonical-form · notation · superseded path · declared tag). → L-059, L-060, L-065, L-071, L-077,
  L-109
## 4. Retirement & staleness: the unit of repair is the THESIS

The three searches, the surfaces a symbol grep cannot reach, marker migration, and what a retirement
does to the surviving gates — including the DEMOTION of a rewired comparison and the doc re-scoping
it forces (D.14) — are `retirement-audit` A–G. What follows is what that skill lacks. → L-044

- **Sort every surviving mention by TENSE into FIVE registers, one repair each** (`retirement-audit`
  B.19 gives the first two): live guidance ⟹ re-word to the successor, which is a CHOICE you must
  measure · history ⟹ prose stays past tense, only the role is downgraded to a `` ``literal`` `` ·
  landed-but-written-as-future ⟹ re-tense in place + one dated note, never delete the bullets (the
  costliest and least greppable — it reads as a plan, not a claim) · aspirational-but-refuted ⟹ ⛔
  *closed as NOT APPLICABLE* with the structural reason · an ADDRESS (an anchor or eq-label carrying
  the retired word) ⟹ KEEP, and say why. A SIXTH register is the genuine grey zone — present-tense
  GRAMMAR over a historical SUBJECT — and its repair is neither: **date the spelling in place**
  (*"the two-gain spelling here is Wave O's, which is what the cited captures were taken against;
  `N_{2n}` joined at CS4c step 3 and rides this argument unchanged"*), which keeps the evidence
  honest instead of retro-fitting a current member list onto a measurement that never saw it. And
  the special-case-vs-history fork is CHECKABLE, not a judgement: read the composition site, or
  read the fixture's constructor call, where the ABSENCE of a kwarg is the proof.
  `[M]` one sweep: 14 updated · 32 period-history · 9 address
  · 3 genuine referent. The per-site ladder for a heavily-referenced entry point is the same
  partition. → L-066, L-070, L-072, L-019, L-114
- **A cross-reference inside HISTORY is a category error** — a role claims the symbol exists NOW at
  THAT path, so the surviving CLAIM licenses the role, not the surviving CLASS. Put the rule in the
  page as a head-of-block `.. note::` and corroborate by counting both spellings in the same file.
  One register down: a raw FILE PATH in a literal (40 of 100 catalogue-wide dead) — fix the CLASS by
  stating once that which tests catch an ERR is the marker set, never prose. → L-062, L-061, L-111
- **A FIELD or CLASS SPLIT is not a rename: one survivor inherits the retired name's letter, so
  every pre-split sentence is INVERTED, not merely stale** — and every hit resolves, so no gate can
  rank it. The instrument is an AST census of production CONSTRUCTION SITES per package; the output
  is a SYMBOLS block whose fourth column is *Was*. A split can refute the corpus's own thesis, and
  the repair is sharper than the original. Its sibling: a retirement that REMINTS the freed name
  onto a different live object makes the name a homonym across ONE commit — disposition each mention
  by what the PASSAGE describes, never by the name. → L-091, L-077, L-017
- **When the residue is LARGE and is a SIMPLIFICATION, DECLARE it at the chapter root** (machine
  header + a `.. note::` naming it a deliberate simplification, pointing at the canonical label),
  fix only the sites describing the SHIPPED object, and report the census with its denominator.
  Sweeping 37 pedagogical sites is a numerics adjudication riding inside a docs pass; silence leaves
  37 false sites. → L-077
- **Symbol collisions: inventory every MEANING before renaming, and rename only the meaning with NO
  constituency.** Importing algebra from code imports its collisions — KEEP the code's spelling
  (internal consistency outranks local awkwardness) and pay with a `.. note::` naming each overload
  and its disambiguator. A short suffix collides as badly as a one-letter symbol; a NEW page
  assembled from several sources is the prime site for a within-document collision no build sees.
  → L-051, L-059, L-065, L-011, L-025, L-034
- **Preserve the WHY; tombstone, don't delete.** Retitle to the CONCEPT and KEEP the anchor (a bare
  `:ref:` renders the new title, so every citer improves for free); rename an anchor only when its
  name encodes a REFUTED concept, counting citers first and moving them in the SAME edit. A section
  corrected twice gets a SECOND dated note beside the first. A tombstone is prose you are authoring
  NOW: it owes evidence, and may only assert what YOUR page controls. ⭐ **A deeply-nested phase
  chain has TWO reading paths and both need one**: a reader scanning the block's top-of-section
  step summary meets the stale terminal claim ~700 lines before the deep tombstone renders, so
  place a forward-pointer `.. note::` at the TOP of the phase block listing which terminal
  decisions were reverted, BESIDE the inline tombstones. One close-out section does not catch a
  scanning reader. → L-007, L-015, L-040, L-063, L-076, L-066, L-056, L-103, L-114
- **A scope note is load-bearing only BESIDE the claim it scopes** — a correct explanation one
  paragraph away is read by nobody who lands on the equation (twice a page already carried the
  right note 2–6 lines below the site and it did not count). Moving it up is a genuine
  improvement, not gate-gaming, and it is why an annotation window is specified in LINES. → L-114
- **Discharging a seam edits FOUR surfaces, and the one nobody edits is the section's own CARDINAL
  NUMBER** ("**Three** arms are deliberately not built"). Grep the section for its cardinal and for
  every forward-looking verb, not only for the seam's noun; past-tense the WHY in place. → L-075
### 4a. The changelog contract

- **Read the page's own preamble before deferring to what N neighbours do — its exception clause
  decides your case.** `history.rst` contracts *"a new entry lands with its merge hash or not at
  all"*, excepting only a Where naming an UNMERGED BRANCH; `spaces.rst` / `field_algebra.rst` /
  `operator_algebra.rst` carry the `*(in development)*` hatch. Route an unmerged entry to a page
  that permits the hatch and report the blocked row ready-to-paste; never fake a hash. → L-099,
  L-067, L-063, L-057
- **The hatch is a DEBT the merge calls in: on discharge grep the BRANCH NAME across `docs/`
  first**,
  then `in development)\*` corpus-wide — nothing points at those cells, because the dispatch names
  only the blocked page. → L-068
- **A dated row keeps the spelling current on its date**, so its discriminator is *dead xref?
  present tense?*, never *stale spelling?* — a uniform sweep would destroy the as-of-its-date
  record. But a CONSEQUENCE clause inside a dated row can be repealed by a later row: tombstone in
  place naming which half survives, and in a reverse-chronological table *"the row above"* is a
  correct pointer. → L-111, L-106, L-068
- **A changelog's staleness markers ROT, and the cheap tell is a row contradicting itself across two
  columns.** Check the page's chronological DIRECTION before placing (per-page, two lines); a date
  in a prose history block is a git question and drifts by one day. → L-107, L-076, L-068
- **Group a big merge's rows by THESIS, not by the plan's phase labels** — the page's own precedent
  settles one-row-vs-many in one grep, and plan-internal tokens are stripped on the way in. → L-068

---

## 5. Page surgery: slice programmatically, assert before writing

- **Never hand-retype a block.** Read → slice → `"".join` → write, with ALL structural asserts run
  on the in-memory result BEFORE any write, so a failed assert leaves the tree untouched (no
  `git checkout` recovery, which `process-discipline` forbids anyway). Strongest guard for a
  confined edit: BOUNDARY BYTE-IDENTITY (`src[:i] == out[:k]`, `src[j:] == out[m:]`) plus an exact
  length delta; for a bulk delimiter edit, `src.replace('`','') == new.replace('`','')`. ⚠ A red
  guard may be the GUARD's error — diagnose whose failure it is first. → L-012, L-022, L-023,
  L-026, L-058, L-062, L-061
- **PROBE docutils, never reason about it** — a stub-directive harness re-checks a 5 800-line file
  in under a second and settles three questions per call (emphasis ⊃ strong warns; emphasis ⊃
  literal is silent and renders raw backticks; `key=``x``` warns, `key=\ ``x``` is clean). → L-061,
  L-056
- **A uniqueness guard over labels or titles compares EXACT LINES, never substrings** — eq-label
  families are BUILT by suffixing, so prefix collision is the normal case here. → L-060
- **A directive with a rendering body needs a placement rule, or 50 land mid-sentence:** after the
  `.. math::` block, unless the next paragraph is a grammatical continuation (`where …`, `so …`,
  `with …`). Read the spanning sentence out loud after writing any such directive. → L-060, L-071
- **⛔ A block replacement that ends MID-PARAGRAPH is swallowed by the directive you insert** —
  extend `old` through the sentence, or re-emit the tail at body level. → L-101
- **Locate by STABLE TITLE, never by the brief's line numbers**, and prove contiguity by counting
  ALL H1 underlines in the range (an anchorless sibling H1 is invisible to the anchor-grep the
  brief's author used). Splice mechanics that broke builds: a slice joined directly before a
  `.. _anchor:` GLUES it (join with `\n\n`); re-nesting demotes every migrated underline
  length-preservingly; removing a middle H1 auto-reparents its H2. → L-022, L-026
- **A cross-page MOVE is ref-safe if the labels MOVE with it** (defined exactly once); only `:doc:`
  needs fixing, plus now-intra-page `(:doc:sibling)` parentheticals, which lie without warning.
  Labels are PATH-IMMUNE — the real break is consuming prose naming the old page, swept with a
  whitespace-FLATTENED scan. → L-022, L-024, L-026
- **Metadata relocation, not deletion:** strip campaign provenance into the changelog, KEEPING
  invariants, eq-labels, vv-status, gotchas, `#N` refs and numerical data, mapping each item to a
  destination FIRST. An overloaded-symbol sweep re-classifies EVERY survivor in a final grep.
  → L-028, L-011, L-025, L-034

---

## 6. Match the doc SHAPE to the event class

One line per event class; the shape itself lives in the archive section named.

- **A rename ruled on a DOC's own argument owes that doc a tombstone** — keep the ruling word for
  word, re-point its subject, spell the retired name as a literal. → L-112
- **An ERR entry is earned by a defect of SILENCE, not of ABSENCE** — the discriminator is *can it
  be re-introduced and reddened?* Its best content is the HIDING mechanism (a COUNT gate over an
  inventory pins whatever inventory it was written against, so an omission reads as intentional).
  → L-110
- **A doc's ARITY is an API.** Before adding a tier to an N-tier table, read the table's own arity
  claim and its guard note, and grep the corpus for the numbering (37 lines cited `Layer-1…4`,
  including an unrelated homonym). Prefer the page's own sub-lettering to a renumber, and give
  genuinely new objects their own H2 saying *these are not a fifth layer*. → L-109, L-104
- **A COUNT-WORD names a SET, and adjacent pages mean different sets by one numeral** ("five
  operators" `{L,C,S,B,F}` beside "four operators" `{L,C,S,B}`). Write the MEMBERS beside the
  count, or rename to a countless form; then triage every count-word hit BY REFERENT — of 29
  `(four|five|six)[- ](term|operator)` hits, most were a different four. ⭐ The nastiest member is
  a number matching nothing on the page: **do not bump it — re-derive what it counts**, and drop
  it if nothing does. → L-114, L-109
- **When a pass changes an OPERATOR, grep its SPLITTING, its ITERATION MATRIX and its
  PRECONDITIONER spelling too** — `ψ_{n+1} = (L+C)^{-1}(Sψ + Bψ + q)` carries no `L+C−S−B`
  substring, so a spelling census is blind to it, and leaving it makes the page state a
  five-member operator whose own splitting drops a term. Two more classes the census cannot see:
  a SECTION HEADING naming the count, and a non-rendered machine-facing
  `.. (vv-status rationale)` COMMENT restating the retired member list. → L-114
- **A LEDGER gaining a field splits across two pages by REGISTER** — theorem and admission grid on
  the point-set page; ledger, SYMBOLS block, per-geometry derivation and worked examples on the
  algorithm page; each cites the other once. → L-091
- **A KERNEL CHANGING HOUSE splits into two registers** — the mathematics stays where its argument
  already flows; the module's own section owns the BOUNDARY (why the verbs are the measure's, why no
  façade, the call-site proof of one closure, the numerical evidence). "ONE closure" publishes as a
  CALL-SITE COUNT; an inert architectural step is named inert WITH its denominator. → L-090
- **A BRANCH-BECOMES-ONE-FORMULA carve has five moves** — one labelled equation in one home and
  every other site points · retitle the section whose title states the refuted claim, keeping the
  anchor · a `.. note::` saying precisely WHAT WAS LOST · a name with three generations as one
  bullet list · the Mode-12 blindness as a labelled subsection with numbered consequences. → L-089
- **A NEW LAYER nobody owned gets its OWN page**, and the decisive argument is
  self-undermining-if-homed-elsewhere. Manage the twin risk actively: own only the register that
  measures 0 hits corpus-wide, opening *"Edited there, consumed here."* → L-079
- **A MINT WITH ZERO CONSUMERS states ⛔ *capability, not a fix* three times** — Key Facts, the seam
  table's first row, and every page touched, each with its own `[M]`. → L-079
- **A ROUTE RE-POINT earns its own labelled doctrine section**, because every value gate over it is
  `X == X`; publish the observation, the DECOY instrument, and the per-decoy admissibility table
  with the MECHANISM behind each measured floor. → L-075
- **A MODELLING-TRUNCATION correction is ONE anchored `.. warning::` carrying the whole measurement
  set; everything else POINTS** (physics home · data-layer home · Key Facts · machine header, one
  clause each). An anchor above an admonition needs explicit link text at every citer. → L-078
- **A REFUSAL BECOMING A CAPABILITY has five moves, and move 4 stops the over-read** — keep the
  diagnosis and past-tense only the verdict · the refusal era verbatim under a ⛔ with why it was
  correct then · SPLIT the recorded debt · say what did NOT ride along, measured · a new ERR
  CHAPTER, not a new number, since the landed gates already carry the old `catches`. → L-076, L-065
- **A NEW THEOREM is homed where its LOCAL half is already derived; then audit the UNIVERSAL it
  amends tree-wide**, where the adjudication is NOT uniform (scope-to-bulk · tombstone · one clause
  · LEAVE). Measure every gate fixture the corpus names against the new predicate. → L-057
- **A POSING-CONTRACT section has six parts in order** — the fields and why there is no default ·
  why the substrate kept what it kept · the guard ruling as an attack table · the ROUTE gate · the
  performance ruling verbatim with a cost table · *What moved, concretely*. → L-072
- **An UN-WELD doc's load-bearing content is the FORCING, not the twin** — quote the layer contract
  that made the duplication unavoidable, and pair it with the honest scope, since the headline is
  nearly always over-broad. → L-071
- **A STRUCTURAL "can never apply" claim publishes as an IFF with numbered conditions + a per-family
  adjudication table + a `.. note::` "what WOULD change this answer"** — the last is what pre-empts
  the "so it's just not built yet" re-reading. → L-070
- **A DIALECTICAL SEED PAGE is not the 9-step close-out arc** — Key Facts carrying the discriminator
  tests verbatim → taxonomy → theorem → the doctrine dialectically (each refutation titled with the
  REFUTING QUESTION, not the verdict) → fences → dev history; and say what the doctrine does to the
  tension it settled. → L-064
- **A fix that works BY RETIRING a failed-approach family gets a SUCCESS-RESOLUTION chapter, not the
  9-step CLOSED arc** — one supersession banner plus targeted tombstones on bald reversals only;
  flip any prior close-out's "open research path" that landed. → L-013
- **An ONTOLOGY-OVERTURN changelog goes on the page whose THESIS moved**; on a blocked page,
  tombstone only the falsified HALF of its row. → L-063
- **A correction sweep must not acquire a SECOND SUBJECT** — fix the claim you were sent for, REPORT
  the neighbour with its proofs. The boundary: co-fix a neighbouring defect only where it sits
  INSIDE a clause you are already rewriting; a standalone instance of the same family is flagged
  for its own pass, with its sites listed. → L-078, L-114
- **A page that says "the OTHER method does X instead" is the sentence nobody re-reads** — it lives
  on page A and its truth lives in solver B, so neither method's maintainer revisits it (one such
  had been retired by an ERR entry two months earlier and was still present tense). ⟹ open the
  other solver's function AND grep the error catalogue for the term; the catalogue entry is where
  a reversal is recorded. Repair: a `.. note::` keeping the retired sentence verbatim in quotes,
  naming what retired it, and stating the surviving difference. → L-114
- **An ARCHITECTURE pass has no Branch-1 SymPy, and its source ladder is the OTHER one:** rich CODE
  docstrings (module / class / property) → the throwaway instruments in `derivations/diagnostics/`,
  which you READ *and RUN* to confirm every cited number → the commit bodies, which carry the exact
  ULP figures and the rationale, QUOTED → the cross-domain-attacker frame memos, where the hardest
  "why" lives. These ARE the algebra of record there. Pull numerical bounds from the TEST FILES,
  never from a brief's memo estimate. → L-114
- **A campaign landing step after step into ONE page uses a shared 8-part peer `====` template:**
  lead (commit chain · issue · date · one line placing the step against its predecessor) · Key
  Facts, 5–7 bullets, always one honest-scope · labelled equations each with its `vv-status` · ONE
  list-table that is the canonical at-a-glance index · the load-bearing rationale, usually a
  rejected-alternative catalogue (the Cardinal-Rule-3 payload) · a numerical-evidence list-table ·
  an honest-scope `.. warning::` framed as the attacker's ABSENCE, not a hedge · cross-refs to the
  predecessors. The per-solver page gets Key-Facts bullets pointing here, never a second
  derivation. → L-114
- **Other event classes, one pointer each:** campaign capstone arc → L-039 · deepening a documented
  feature, a PLANNED-not-built admonition paired with a current-state subsection → L-014 · an
  EVICTION changes the CARRIER, not the physics → L-016 · a completed architecture earns ONE
  taxonomy-culminating section → L-018 · a NEW foundational chapter reads then RUNS the
  algebra-of-record module first → L-025 · growing a thin honest stub at campaign close → L-036 ·
  "is the terminal docs phase done?" almost always answers effectively-DONE, and a documented SEAM
  is the opposite of a gap → L-038 · merging a re-staged branch's docs into a diverged tree → L-012.

---

## 7. V&V vocabulary — you are the curator (AGENT.md Directive 5)

You write the prose future readers QUOTE about verification status. The level definitions, the three
pillars and their evidence boundaries, the anti-patterns and the failure modes are `vv-principles` —
match it VERBATIM, never paraphrase. → L-010

- **When the WHOLE point is that the obvious reading is WRONG, name the FORBIDDEN SENTENCE** in an
  explicit `.. warning::`, verbatim, with the correct framing beside it. A future session quoting
  the page for V&V reasoning meets the warning before the misreading. It is the page-level
  analogue of `vv-principles`' *"NEVER write 'MMS verifies the eigenvalue'"*, and it beats a hedge
  because it is greppable. → L-114
## 8. Code-prose rebalance (docstring / comment trimming)

- **Expect ZERO MOVED.** Cardinal Rule 3 means the theory shipped WITH the code, so a concept that
  feels unique to a file is almost always already TWIN in the landing chapter — grep the chapter
  before crediting one MOVED; a pre-classifier's MOVED column is ~100 % noise. → L-033
- **The CONTRACT test: "would a competent modifier who never leaves this file do the wrong thing
  without this line?"** If yes it is CONTRACT, however history-flavoured. → L-033
- **FILE-CLASS sets the size and the SURFACE of the honest cut** — teaching-heavy operator ⟹
  aggressive twin-cut; contract-heavy operator / machinery / ABC ⟹ small, and the surface is
  module-head essays and duplicated numbers, not method bodies; driver / mesh ⟹ hunt standalone
  `#`-comment tombstones first. A −2 to −5 % cut is CORRECT: report the file-class rationale so it
  is not read as timidity. → L-034
- **Provenance trimming is citation-vs-narration applied uniformly** (trim landed campaign-STEP
  codes, KEEP bare `#NNN` anchors and named patterns) — but a hand-transposed-adjoint comment body
  IS the algebra of record. A batch "special" is a VERIFICATION obligation first, an edit obligation
  only on failure. → L-034
- **Prove the edit is doc-only by AST/token comparison against HEAD**, not by reading the diff — it
  also proves no `verifies`/`catches` marker moved. It is blind to comments (fine) and an f-string
  assertion message is CODE: leave it and REPORT. Run the Sphinx gate iff the file is
  `automodule`'d; `:noindex:` does not exempt it. → L-041, L-045, L-073, L-033

---

## 9. Gates, generated artefacts, tooling

- **Generated artefacts are NEVER hand-edited** (the V&V matrix, capability tables,
  `_generated/*.inc.rst`, `vv-principles/error_index.md`) — fix the registry-side metadata and
  report
  the REAL post-regen number. A `-E` build on a dirty branch absorbs rows from other uncommitted
  work: report it, never revert it. In a fresh worktree a missing generated artefact is an ENV gap,
  not a doc defect. Orphaned built HTML from a renamed source looks like a live stale ref —
  discriminate by "does the source `.rst` still exist?". → L-008, L-026, L-040
- **⛔⛔ EVERY corpus grep owes `| grep -v _build`** — four generations of stale HTML answer for the
  source: a briefed anchor returned 20 hits, ALL under `docs/_build/`, where the source count was 0.
  A brief saying "add X where Y is" owes a check that Y is there. → L-104
- **⛔⛔ An `.. error-entry::` without its `catches` marker REDDENS the suite** — a docs-only pass can
  break it, there is no machine-readable exemption, and the generated index moves too. Ship the
  entry
  ready-to-paste with both marker lines and the build order, and report the marker as BLOCKING.
  → L-091, L-103
- **ONE re-runnable python self-check beats the build for structure, at ~2 s:** short-underline
  detection, ladder-order, per-table column consistency, `:widths:` sums, label/anchor uniqueness by
  exact-line compare, role import-resolution, and `:eq:`/`:ref:`/`:doc:` resolution against the
  whole
  corpus. It caught every structural defect on a 1 158-line new page before any build ran. → L-064
- **Measure a gate baseline from `git archive HEAD` into a temp tree** for a true before/after on a
  dirty working tree — and read its traceback before counting its warnings: untracked DATA files are
  absent from the archive, so the pristine tree can carry artefacts the live tree does not. `rm -rf`
  inside a compound Bash command is refused here; `mkdir -p <fresh dir>`. → L-059, L-073, L-051
- **An error-message string inside `raise` is EXECUTABLE — report it, don't edit it** under a
  doc-only constraint. ⭐ And check WHICH SUBSTRING is pinned first: tests match the shortest
  distinctive fragment, so "pinned, leave it" and "unpinned, safe to correct" are different
  instructions to the code owner (`retirement-audit` B.11). → L-041, L-073
- **A corpus-wide mechanical migration is dry-run-first and WHITELIST-scoped**; key any block
  remover
  to INDENTATION too, or it eats footnotes. → L-031
- **Self-check the V&V scan directly, not via the full audit** — `_scan_theory_equations` runs in
  <1 s, avoids pytest collection, and does not trip on a sibling batch's in-progress sentinels.
  → L-035, L-063, L-069
- **Widening someone else's issue: re-run THEIR instrument, not yours** — their number plus your
  wider denominator plus the exclusions retitles the issue; a fresh regex forks the count. → L-068

---

## Quality self-assessment (AGENT.md Directive 3)

Rate the six dimensions and log the weakest. On TERMINOLOGY / ROUTING / retirement passes the weak
dimension is routinely "numerical evidence" — structurally ABSENT (no flux moves ⟹ no convergence
table), not a deficit. Say so; do not manufacture one.
