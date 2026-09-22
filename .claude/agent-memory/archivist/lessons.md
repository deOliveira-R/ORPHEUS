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

- **A relayed claim is a PREDICATE, not a fact — read the slot in the live tree before re-wording
  it.** "renamed everywhere" licensed the inference that every slot holding the hub was renamed;
  24 old-spelled params survived, and "fixing" the quoted live signatures would have read as the
  careful half. Publish the residue with its count and predicate. → L-112, L-001, L-018, L-039
- **A briefed LIST can contradict the brief's own carve-out; that contradiction IS the finding.**
  CONSTRUCT each member and print the property — a grep finds the members that have it and is
  silent on the ones that do not, which is the half that matters. → L-096, L-095, L-104
- **A brief's named target can measure ZERO — run the FILE's own predicate before reporting
  "nothing to do".** A 0 on someone else's grep is not a clean file, and a phase that landed with
  no docs pass leaves its rot for the next phase's sweep: budget for it. → L-102, L-072, L-075,
  L-081, L-104
- **A brief can instruct a HARD audit error — read the scanner before obeying a `vv-status`
  instruction.** `documented` is the ONLY legal status (`tests/_harness/audit.py`, exit 2
  otherwise); what a new `verifies` marker earns is UN-SENTINELING, not an upgrade. → L-106
- **A brief describing code being typed NOW is a FORECAST** — verify the symbol, and poll on the
  INVARIANT (`! git grep -q <old> -- orpheus/`), never on one new symbol: a half-written tree
  answers plausibly, and its `AttributeError` reads as "unsupported", not "broken". → L-087, L-098
- **A brief's discriminator is a heuristic; resolve EACH site by its own live signature.** A
  mechanical vocabulary swap restates a false claim in fresh, authoritative words. → L-001, L-021
- **A brief's "sharpest observation" is a HYPOTHESIS with a computable confusion matrix** — publish
  the measured split and the refutation, or the next reader re-derives the heuristic. → L-060
- **A design review's verdict is a RECOMMENDATION, not what shipped** — `dataclasses.fields` the
  class before quoting its field list, and publish the divergence as a refuted-candidate row with
  its reason (`process-discipline` "A refuted candidate is first-class output"). → L-104, L-001

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

- **Re-derive every numeric literal the pass publishes, in ONE script, at the end.** A spot check
  finds neither of the two defects the literal SET finds (a figure I had invented; a figure that was
  real and answered a different question). → L-109, L-001, L-042 · now AGENT.md checklist item 6
  (2026-09-21)
- **Two honest `[M]`s of "the same thing" differ by the PREDICATE or the STATISTIC — never
  adjudicate; find the reading that makes both true and PRINT the arithmetic.** An exact
  reconciliation proves the census complete AND names what it excluded; when nothing reconciles,
  RETIRE the pair and re-measure with the fixture in the caption. → L-107, L-108, L-100, L-093,
  L-085, L-091, L-051
- **When a quoted figure will not reproduce, find the PARAMETER and publish the CURVE** — a single
  figure without its tolerances is one point on a curve wearing a constant's authority. → L-110,
  L-108, L-065, L-067
- **Publish the MECHANISM or the CLOSURE ARGUMENT, not the value.** A derivation stays true as the
  tree grows; a number rots and gets COPIED into docstrings. Instances: a mutation magnitude is a
  draw (publish the absorbed piece, not the digits) · a published digest LENGTH is an artefact of
  your throwaway class's name (publish *byte-EQUAL*) · a positive control's VALUE carries no part of
  the argument (publish the PREDICATE, keep the zeros). → L-064, L-104, L-103, L-070, L-075, L-088
- **A `[M]` whose denominator is a COMPUTED SET has a shelf life the FINDING does not** — say "the
  finding is unchanged and only the DENOMINATOR moved, because it is the size of a candidate set".
  → L-090, L-091
- **Don't assert a mechanism — MEASURE it.** "A flat vector broadcasts across the spatial axis" was
  true and useless; measured, it is an outer product, which buys three traps the hand-wave misses.
  → L-109
- Census and number discipline is `instrument-doctrine` X2 with `plan-authoring` §2/§4. The three
  clauses this corpus breaks most: **QUANTIFIER** (every "each/every/all" you publish is countable
  in one command, and the measured sentence was strictly better every time) · **DRAW**
  (bit-exactness, a ULP gap and a percentage are properties of the draw — publish a bound over
  ≥200 seeds and say whether you have construction-exact or draw-exact) · **RELAY** (publish YOUR
  number with YOUR configuration). → L-060, L-067, L-071, L-076, L-084, L-099, L-057, L-043

### 1d. Read the object, not its description

- **Construct the object and print the property.** `hasattr` is False for a dataclass field with no
  class default and for an `__init__`-assigned attribute; the ladder is `hasattr` →
  `dataclasses.fields`/`__annotations__` across `__mro__` → `self.x =` across the MRO → CONSTRUCT.
  → L-093, L-076, L-053, L-105
- **For "all N of family F do X", re-derive F's MEMBERS, not |F|** — a family turns over with its
  count fixed, and a roster's own list is not the family (`vars(cls)` +
  `isinstance(v, classmethod)` said FIVE where four pages said four; the fifth row was the
  argument). → L-101, L-074, L-081
- **A composite's measured identity cannot certify its FACTORS**, and an operator-movement claim has
  a draw-free form: build the matrix column by column with `e_k`, never one probe vector. → L-049,
  L-076, L-084
- **The pre-carve tree is measurable** (`git worktree add … HEAD --detach`, or `git archive HEAD`),
  and when the retired object is shape-characterisable a hand-built STAND-IN beats a diff: the
  equal-shape CONTROL is the finding, and a diff cannot produce it. ⚠ the venv's editable install
  hooks `sys.meta_path` and OUTRANKS `PYTHONPATH` — strip the finder and print `orpheus.__file__`
  as proof, or use a RENAMED shadow package. → L-095, L-101, L-050, L-083
- **Any cross-carve table owes a line naming WHICH statistic is comparable** — the tell that you
  need it is that the CONTROL column moved. → L-095

### 1e. Surfaces nothing gates

- **A published `.. code-block:: python` is the highest-severity staleness there is** — it promises
  reproducibility and nothing gates it. After any constructor/signature change grep the code-block
  BODIES before the prose sweep, and RUN what a present-tense sentence promises works. → L-077,
  L-047, L-093, L-099, L-105
- **A docstring quoting another docstring, and a doc quoting a production docstring, are ungated
  cross-file dependencies — grep the quoted FRAGMENT, not the symbol.** The quotation goes false
  while every name in it stays alive, and it is usually load-bearing for the neighbour's
  justification. → L-102, L-066, L-062
- **A brand-new docstring can ship a convention the code does not have** — it documents the DESIGN.
  Document the VALUE, date it, report the mismatch; never promise the planned flip. → L-109
- **A table's CAPTION owns its columns**, and in a historical table exactly ONE column is
  present-tense — the only one that can rot. "Update column X" can be an instruction to falsify the
  caption; place by CONTENT, never by a site ORDINAL, and write the deviation into the report.
  → L-098
- **A regenerated cache or untracked artefact is a doc surface with no symbol in it; `ls -l` is the
  only instrument** — `git status`, `-W`, the xref gate and `dead_references` are all blind. → L-093
- **A page can contradict itself 80–200 lines apart, and the stale half is the one a reader
  quotes.** After an algebra change grep the OLD spelling WITHIN each page that already carries the
  new one; when a page contradicts itself, the HEDGE is usually the true half. → L-077, L-078, L-092
- **A page's own ⚠ CAVEAT can BE the diagnosis a later carve lands, or the theorem the next step is
  built on — keep it verbatim under a ⭐, never tense-flip it.** Tell: *"X is not always Y, so the
  precise statement is Z, not a named attribute"* is a page saying where its subject is welded.
  → L-104, L-086

### 1f. Landings, predictions and seams

- **A landing is audited part by part, never by flipping a tense.** A deferral row has a CLAIM, a
  PRECONDITION and a MECHANISM; a prediction names a MECHANISM, a HOST/TYPE, a PHASE and a
  DELIVERABLE. Any subset can move: say which part survived, annotated under the row. → L-107,
  L-100, L-049, L-082, L-087, L-097, L-110
- **Keep a refuted prediction — the refuted MECHANISM is the interesting half.** "The defect was
  correctly identified and the mechanism was not" is the durable sentence; a tense flip destroys it.
  → L-108, L-090, L-074
- **A two-clause ruling may have an occupant for only one clause** — publish clause 1 as the ruling
  and clause 2 with the census showing it empty. → L-108
- **A campaign day lands SIBLING steps, and each repeals a PREMISE elsewhere — grep the PREMISE and
  the EVENT-phrases** (*"when the hub gains"*, *"until X lands"*, *"deferred to"*, *"becomes real
  when"*), never the step's name: every symbol is alive and the sentence is grammatical, so no
  build, `dead_references` or symbol grep sees it. → L-081, L-107, L-075
- **When a landing gives a REAL NAME to a phrase the page used loosely, the phrase becomes a
  MIS-STATEMENT, not a head start** — grep the new field's name across `docs/` and read every hit as
  a claim about the new thing even where it predates it. → L-074
- **A landed change that SILENTLY PRESERVES a published table is itself publishable** — re-measure
  and add a dated note saying it survived BY DESIGN, or the next reader assumes staleness. → L-074
- **A `==`-gate's blindness is a property of the IDENTITY RELATION, not of the gate** — a re-typing
  upgrades or downgrades every equality assertion in the corpus with no file touched. Afterwards
  grep `metric-blind` / `cannot tell them apart` / `(name, shape)`-equal and re-derive each: some go
  false, some become STRONGER claims that now under-sell a real gate. → L-100, L-103

### 1g. Instruments and corroboration

- **Your own gates fail silently and flatteringly** (`instrument-doctrine` X1). Two corpus
  corollaries: a gate that can print 0 on an empty input list must echo the INPUT COUNT (an `instrument-doctrine` skill X1 clause since 2026-09-21); and a
  LINE-based scanner cannot see a role, a bold run or a claim split across two lines — scan
  CONTIGUOUS blocks with `re.S`, and use RETIRED symbols as positive controls. → L-108, L-107,
  L-101, L-088
- **Two prose surfaces agreeing is ONE surface** (X4), and git is the arbiter for dates and merge
  status (`process-discipline` "Trust git for merge status"). Two surfaces agreeing on a NUMBER are
  one surface until you check they do not share a seed. → L-064, L-099
- **Two independently-VOCABULARIED instruments agreeing IS the acceptance evidence** — and adopting
  a reviewer's exact patterns is evidence, not concession. → L-052, L-067, L-070, L-086
- **You are the judgment layer over any bulk scanner** — import-verify every suggested target,
  reject findings the clean build contradicts, attribute every residual to a named false-positive
  class. → L-021

---

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
- **⛔⛔ TWO STANDING CLAIMS OF MINE ARE REFUTED, both by one 20-second throwaway project with
  two-sided controls** (`[M]` 2026-09-21, Sphinx in this venv, a 2-page project, live targets beside
  dead ones). **(a) A cross-doc dangling `:ref:` DOES warn** — `WARNING: undefined label: '<label>'
  [ref.ref]` at DEFAULT severity, bare AND with explicit text, exactly like the intra-doc case; a
  dangling `:doc:` warns as `[ref.doc]`. So "a cross-doc `:ref:` miss is silent at every severity"
  is FALSE, and every label-rename caution argued from it is void: renaming a section anchor is
  build-gated for its `:ref:` citers just as an eq-label is for its `:eq:` citers (L-070). **(b)
  `-n` DOES catch dead python-domain roles ON AN `.rst` PAGE** — all five of
  `:class:`/`:func:`/`:meth:`/`:attr:`/`:mod:` warned under `-n` and NONE at default severity. The
  true scope of "`-n` does not save you" is the surface Sphinx never RENDERS — a docstring in an
  un-`automodule`'d module, and everything under `tests/` (L-044's measurement, correctly). ⟹ for
  `.rst` pages `-n` is an available instrument for the exact class the project xref gate is blind
  to (§2b; `retirement-audit` A.2 carries it since 2026-09-21); use it as a pre-edit-vs-post-edit SET DIFF, because plain-text-by-convention roles into
  the 24 noindex modules nitpick as "target not found" by design. → L-113, L-002, L-044, L-063, L-070,
  L-076
- **⛔⛔ `:noindex:` on an `automodule` mints NO cross-reference target, so the role renders plain
  text.** `[M]` 2026-09-21 (mine, by parsing each directive's option block): **48** automodules in
  the source, **24** of them `:noindex:`. `[M]` 2026-09-18 (L-112, not re-run): a two-sided HTML
  control gives 0 `id=` anchors for a noindex module beside a normal page-mate, and **1 081**
  python-domain roles point into the noindex set — re-measure that one before quoting it. REPORT the
  finding; flipping `:noindex:` risks duplicate-object warnings and is an architectural decision.
  → L-112, L-053, L-002
- **These DO warn, so they are gated:** a dangling `:eq:` (measured with controls) and — per the
  refutation above — a dangling `:ref:` and `:doc:`, cross-doc included, so BOTH label classes are
  rename-gated by the build; a bad `:by:`; an unresolved
  `catches("ERR-NNN")`; `ref.ref` *"A title or caption not found"* (a bare `:ref:` to an anchor
  above a paragraph, an admonition, or a **bold run-in heading**); a `:widths:`/column mismatch; a
  malformed `===` table (never hand-align one holding a `:math:` role — source length is not
  rendered length); an italic run interrupted by a role. ⚠ `SyntaxWarning` and the nexus per-marker
  lines carry no `WARNING:` prefix, so read the whole log and grep `CRITICAL:` too. → L-070, L-060,
  L-002, L-027, L-040, L-054, L-055, L-048, L-095

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
- **Plain-text roles are often the page CONVENTION, not a defect** — match the page, repoint dead
  refs to the LIVE path, never half-surface one or two leaves of an un-surfaced package.
  `automodule` readiness is multi-gate ("0 `:label:`" is necessary, not sufficient): `-E -W` each in
  isolation and report the unblocking fix. → L-002, L-034, L-047

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
  a CONTEXT SHIFT, not two events; a page-wide count indicts someone else's prose. ⚠ `nohup … &`
  inside a background Bash call reports the SHELL's exit, not sphinx's. → L-072, L-088, L-089
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
- **Section-label and equation-label are DIFFERENT namespaces** coexisting under one name with no
  warning; verify with `grep -c '^\.\. _X:'` / `grep -c ':label: X'`, never a mention count. The
  WHERE-list under an equation is the tell that the equation drifted from its own prose. → L-024,
  L-003, L-056

---

## 4. Retirement & staleness: the unit of repair is the THESIS

The three searches, the surfaces a symbol grep cannot reach, marker migration, and what a retirement
does to the surviving gates — including the DEMOTION of a rewired comparison and the doc re-scoping
it forces (D.14) — are `retirement-audit` A–G. What follows is what that skill lacks. → L-044

- **The load-bearing half of an ontology overturn has NO dead symbol in it.** Grep the retired
  symbol to FIND the sites, then read the enclosing ARGUMENT to decide the edit: a five-obstruction
  proof can rest on a premise that is now false while its CONCLUSION still holds. Re-derive from
  what survives, keep the conclusion, tombstone the example — the live tree usually hands you the
  replacement. When the deletion is a COROLLARY of a design unification the SECTION'S THESIS is
  stale, not the line; the tell is a stale design stated in the present tense as live rationale.
  → L-063, L-069, L-020, L-013
- **Sort every surviving mention by TENSE into FIVE registers, one repair each** (`retirement-audit`
  B.19 gives the first two): live guidance ⟹ re-word to the successor, which is a CHOICE you must
  measure · history ⟹ prose stays past tense, only the role is downgraded to a `` ``literal`` `` ·
  landed-but-written-as-future ⟹ re-tense in place + one dated note, never delete the bullets (the
  costliest and least greppable — it reads as a plan, not a claim) · aspirational-but-refuted ⟹ ⛔
  *closed as NOT APPLICABLE* with the structural reason · an ADDRESS (an anchor or eq-label carrying
  the retired word) ⟹ KEEP, and say why. `[M]` one sweep: 14 updated · 32 period-history · 9 address
  · 3 genuine referent. The per-site ladder for a heavily-referenced entry point is the same
  partition. → L-066, L-070, L-072, L-019
- **A cross-reference inside HISTORY is a category error** — a role claims the symbol exists NOW at
  THAT path, so the surviving CLAIM licenses the role, not the surviving CLASS. Put the rule in the
  page as a head-of-block `.. note::` and corroborate by counting both spellings in the same file.
  One register down: a raw FILE PATH in a literal (40 of 100 catalogue-wide dead) — fix the CLASS by
  stating once that which tests catch an ERR is the marker set, never prose. → L-062, L-061, L-111
- **A stale REASON outlives a stale NAME, and only the name is greppable** — after fixing a retired
  name, read the sentence that JUSTIFIES it: keep the imperative, replace the *because*, say what
  changed. A retirement also propagates to BOUNDS, to NEGATIVE claims (`independent of|unaffected
  by|does not depend on`) and to WORKED EXAMPLES and EXHIBITS, none of which a symbol grep sees.
  → L-069, L-067, L-054, L-082
- **Fixing HALF a claim in one file is worse than fixing none.** After repairing a section, grep the
  WHOLE FILE for every spelling of the retired predicate and adjudicate each by tense: a
  self-contradicting page is citable for EITHER sentence. A phase-N pass also inherits
  phase-(N−1)'s falsifications — audit the PARAGRAPH FAMILY, not the commit's diff. → L-056, L-042,
  L-077
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
- **Before repairing a stale equation, census the corpus for a page that already states it right.**
  That stops the repair minting a TWIN (keep the equation + an `.. important::` naming the SSOT and
  which claim THIS page owns) and it tells you whether you are fixing an outlier or inventing a
  convention. When the corpus states one object N incompatible ways, each internally consistent, a
  hidden PARAMETER is unnamed — name it ONCE in a table and make every site a pointer. → L-060,
  L-065
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
  NOW: it owes evidence, and may only assert what YOUR page controls. → L-007, L-015, L-040, L-063,
  L-076, L-066, L-056, L-103
- **A CAPABILITY FLIP stales DEFERRAL CONTRACTS, and the blast radius is the WHOLE page.** Census
  with a ±3-line CO-OCCURRENCE window, never a line grep, and publish the PREDICATE, not the count
  (the co-occurrence count RISES when you succeed, because a correction names what it corrects).
  Grep every future-tense token (`blocked|not built|not yet|pending|in flight|future seam|lands
  with`) — a `*(in development)*` hatch WRAPS, so a line-based grep reads 0. The signature defect is
  a STALE HEADER over an ALREADY-CORRECTED BODY: read ±30 lines and adopt the neighbouring truth's
  spelling verbatim. It also stales "X is UNTESTABLE" and "X is unreachable", and when the LAST
  unbuilt X lands, RE-POINT the sentence to the still-unbuilt sibling rather than dropping X.
  → L-073, L-037, L-075
- **Discharging a seam edits FOUR surfaces, and the one nobody edits is the section's own CARDINAL
  NUMBER** ("**Three** arms are deliberately not built"). Grep the section for its cardinal and for
  every forward-looking verb, not only for the seam's noun; past-tense the WHY in place. → L-075
- **A plan-internal STEP LABEL collides like a bare `#N`** — this corpus carries four `R7`s.
  Disambiguate by CAMPAIGN at every use and name the ones you are NOT (`process-discipline`'s
  plan-number clause, at corpus scale). → L-104, L-067

### 4a. The changelog contract

- **Read the page's own preamble before deferring to what N neighbours do — its exception clause
  decides your case.** `history.rst` contracts *"a new entry lands with its merge hash or not at
  all"*, excepting only a Where naming an UNMERGED BRANCH; `spaces.rst` / `field_algebra.rst` /
  `operator_algebra.rst` carry the `*(in development)*` hatch. Route an unmerged entry to a page
  that permits the hatch and report the blocked row ready-to-paste; never fake a hash. → L-099,
  L-067, L-063, L-057
- **`git merge-base --is-ancestor` every unstamped row above yours, INCLUDING your own from this
  morning** — three stale rows in one table, twice. A vanished branch means merged; say whether the
  hash was written verbatim in a plan (zero risk) or DERIVED from `git log --first-parent`. → L-105,
  L-106, L-111, L-100
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
- **When the brief says "relocate to X" and X already carries it, the action is DE-DUPLICATION, not
  relocate+merge** — replace with a `:doc:` pointer, merge nothing, FLAG the inversion. Prefer an
  additive roadmap + `:ref:` to the SSOT over copying a table; a FOLD is a MOVE. → L-027, L-029
- **Run a SELF-CONSISTENCY pass on prose YOU authored before the first build** — new prose that
  DECLARES a rule must obey it; ask of every name *which branch of my stated rule does this take?*
  Author heads and intros as pure literals so no f-string layer mangles LaTeX braces. → L-054,
  L-062, L-026
- **Metadata relocation, not deletion:** strip campaign provenance into the changelog, KEEPING
  invariants, eq-labels, vv-status, gotchas, `#N` refs and numerical data, mapping each item to a
  destination FIRST. An overloaded-symbol sweep re-classifies EVERY survivor in a final grep.
  → L-028, L-011, L-025, L-034

---

## 6. Match the doc SHAPE to the event class

One line per event class; the shape itself lives in the archive section named.

- **A NAMING pass adjudicates SENTENCES, never words** — a hub-SHAPED regex cut ~1 100 raw hits to
  ~135 readable ones. The five keep-classes a blanket replace corrupts: genuinely geometric · a live
  API literal · a RETIRED TIER'S name (history) · a QUOTATION of an archived plan row · a
  cross-method tier only half-affected (add a ⚠ instead of renaming). → L-112
- **A rename ruled on a DOC's own argument owes that doc a tombstone** — keep the ruling word for
  word, re-point its subject, spell the retired name as a literal. → L-112
- **An ERR entry is earned by a defect of SILENCE, not of ABSENCE** — the discriminator is *can it
  be re-introduced and reddened?* Its best content is the HIDING mechanism (a COUNT gate over an
  inventory pins whatever inventory it was written against, so an omission reads as intentional).
  → L-110
- **A TABLE OF DEFERRALS is the highest-rot surface on a theory page** — every row is a prediction:
  audit it whenever the owning campaign lands anything, close with ✅ + date, never by deleting.
  → L-110
- **A doc's ARITY is an API.** Before adding a tier to an N-tier table, read the table's own arity
  claim and its guard note, and grep the corpus for the numbering (37 lines cited `Layer-1…4`,
  including an unrelated homonym). Prefer the page's own sub-lettering to a renumber, and give
  genuinely new objects their own H2 saying *these are not a fifth layer*. → L-109, L-104
- **An "X is arm-asymmetric BY DESIGN" argument that gets dissolved: keep the table, flip the verbs,
  add the price** — principled equivalence, never bit-identity. → L-110
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
- **A RETRODICTION table needs a STATUS column IN the row** (`plan-authoring`
  TABLE-ROWS-ARE-CLAIMS, at corpus scale) and a heading naming the real epistemic claim. → L-064
- **Citing an SSOT: name the REGISTER your page owns** — the same fact in a different register is
  not
  a twin. Derive it a third way as UNLABELLED math, cite the SSOT's label, open with *"Edited
  there, consumed here."* → L-064
- **A fix that works BY RETIRING a failed-approach family gets a SUCCESS-RESOLUTION chapter, not the
  9-step CLOSED arc** — one supersession banner plus targeted tombstones on bald reversals only;
  flip any prior close-out's "open research path" that landed. → L-013
- **An ONTOLOGY-OVERTURN changelog goes on the page whose THESIS moved**; on a blocked page,
  tombstone only the falsified HALF of its row. → L-063
- **A correction sweep must not acquire a SECOND SUBJECT** — fix the claim you were sent for, REPORT
  the neighbour with its proofs. → L-078
- **Stub → rich narrative reads memo → production docstrings → tests → SymPy, in that order** — the
  docstrings are the verbatim prose seed, the memo carries the honest interim scope, and on an
  algebra error you dispatch rather than edit (`algebra-of-record`). → L-005
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

- **A doc sentence "gates X and Y pin claim C" IS a coverage claim** — the prose analogue of
  `vv-principles`' marker rule, and a clause of it since 2026-09-21. Justify it by a MUTATION that reddens X and Y, never by topical
  adjacency, and cite PER FIELD, not per topic. Highest-risk moment is replacing a gate you just
  demoted: the nearest-sounding sibling inherits neither scope. → L-047
- **The SAME gate cited for TWO claims can be right once and wrong once** — on "citation of X is
  false", ask *false for WHICH claim* and grep every occurrence before editing any. → L-047
- **Never upgrade a `@pytest.mark.foundation` gate to an L-level in prose** to make a section sound
  better-verified; read the marks and say "software/structural invariant of a discrete construction,
  not an equation claim". → L-040
- **Prose that summarises a gate outlives the gate's assertion — open the gate and read its
  assertions AND its history note.** A docstring containing *"until <phase> this asserted X"* is a
  gate whose prose summaries elsewhere are presumed stale, and naming which half proves the measure
  is read at all versus which half is the ruling turns a wrong one-liner into design rationale.
  → L-102
- **Write the precise object: the Euclidean transpose `Aᵀ` is not the metric Hilbert adjoint
  `A† = G⁻¹AᵀG`**, whatever a campaign colloquially calls it. → L-010, L-034
- **Skill-uplift duty:** propose the `vv-principles` / `error_catalog` / `algebra-of-record` edit in
  your return whenever you meet a published-prose anti-pattern or evidence-boundary case the skill
  does not capture. The skill grows when you feed it back. → L-010

---

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
- **RE-MEASURE the `-E` baseline every session; never quote a recorded number** (it has drifted
  9 → 1 → 0). Diff the WARNING/ERROR/CRITICAL **set**, not the count. A full `-E` rebuild can exceed
  the 120 s foreground cap. When the red is the carve's OWN, the gate becomes EXIT=0, stated in §0
  of
  the report — "count unchanged" would license shipping 13 errors. → L-027, L-029, L-041, L-094
- **SEQUENCE the session so you build TWICE:** baseline `-E -W` → *all* edits → *all* residual greps
  → xref gate → AST doc-only proof → ONE verification build. Every extra build is bought by an edit
  made after launching, so the self-consistency pass (universals, quotations, denominators,
  superlatives, symbol collisions, aspirational rows) runs to EXHAUSTION *before* the first
  verification build. Re-broken to four builds, then to five. → L-054, L-064, L-081 · now AGENT.md
  checklist item 0 (2026-09-21)
- **ONE re-runnable python self-check beats the build for structure, at ~2 s:** short-underline
  detection, ladder-order, per-table column consistency, `:widths:` sums, label/anchor uniqueness by
  exact-line compare, role import-resolution, and `:eq:`/`:ref:`/`:doc:` resolution against the
  whole
  corpus. It caught every structural defect on a 1 158-line new page before any build ran. → L-064
- **A LENGTH-CHANGING rename breaks section underlines — scan, don't wait for the build**
  (`retirement-audit` F.23): code points, underline at col 0, single repeated marker; fix with the
  FILE's own marker char and assert it. Title markers are file-local, so prefer COPYING a proven
  underline to re-counting code points. → L-112, L-009, L-035
- **Validate your OWN parser against a known-good member before believing its negatives**, and
  record
  this corpus's three standing self-check false positives — a `.. code-block:: rst` example carrying
  a duplicate eq-label, a legal EMPTY list-table cell (`^     -(\s|$)`), and a RELATIVE `:doc:`
  docname — or every run re-litigates them. → L-081, L-089
- **⛔⛔ An UNQUOTED heredoc runs COMMAND SUBSTITUTION on every backtick pair inside it** — a
  `code-search` check since 2026-09-21 (a quoted heredoc to a FILE, paths by env var, `chr(96)`, a
  witness assert per pattern). The founding case: four markup patterns collapsed to "match any
  bold", 120 hits on clean prose; the same collapse the other way prints a clean 0. → L-102, L-061,
  L-030.
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
- **Citations: `grep '^\.\. \[Key\]'` before citing or defining** — resolve cross-doc, never
  redefine, match the page's plain-text convention, and verify the pre-existing duplicate-citation
  count is unchanged rather than zero. → L-006, L-025
- **An agreement `[M]` handed to you is a LADDER unless proven flat** — measure it, publish it, name
  the gate's tolerance, and say a finer row must widen it. → L-054
- **Widening someone else's issue: re-run THEIR instrument, not yours** — their number plus your
  wider denominator plus the exclusions retitles the issue; a fresh regex forks the count. → L-068

---

## Quality self-assessment (AGENT.md Directive 3)

Rate the six dimensions and log the weakest. On TERMINOLOGY / ROUTING / retirement passes the weak
dimension is routinely "numerical evidence" — structurally ABSENT (no flux moves ⟹ no convergence
table), not a deficit. Say so; do not manufacture one.
