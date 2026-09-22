# Rules-uplift candidates

Four. Each names the rule or skill it would join, the clause in that file's own `check:` / `tell:`
form, and its founding case. I propose the wording; the orchestrator edits the SOURCE page under
`docs/development/` (the `.claude/` copies are generated — `tools/harness` rewrites them).

Ranked by what they buy.

---

## U1. `retirement-audit` A.2 — the surface, not the flag, decides whether `-n` is a gate; and `:noindex:` is its own plain-text mechanism

**Joins:** `retirement-audit` skill, item **A.2** ("Text-grep the symbol across code, tests AND
`docs/`, then read the surfaces a symbol grep cannot reach; each surface has one instrument"),
source `docs/development/skills/retirement-audit.md`.

**Why.** A.2 currently says a rendered Python-domain xref is gated by `sphinx -W` *"but Sphinx
nitpicks only what it RENDERS, only ~45 modules are `automodule`'d and nothing under `tests/`
renders, so an unresolved xref elsewhere renders as plain text with no warning at any severity and
`-n` does not save you."* Two corrections, both measured this session:

1. **`-n` DOES save you on an `.rst` page.** `[M]` 2026-09-21, a two-page throwaway project with
   live-beside-dead controls: all five of `:class:`/`:func:`/`:meth:`/`:attr:`/`:mod:` at a dead
   dotted target warn under `-n` (`py:<role> reference target not found … [ref.<role>]`) and **none**
   warn at default severity. The `.rst` corpus is rendered by construction, so the "at any severity"
   clause is true only of the surface Sphinx never renders — a docstring in an un-`automodule`'d
   module, and everything under `tests/`. That is the real scope, and it is the difference between
   "no gate exists" and "a gate exists and nobody runs it".
2. **`:noindex:` is a second, distinct plain-text mechanism, and it is now the majority.** A.2's
   count ("only ~45 modules are `automodule`'d") measures the wrong thing on its own: `[M]`
   2026-09-21 the source carries **48** `automodule` directives, of which **24** carry `:noindex:`
   and therefore mint no cross-reference target at all. Half the "surfaced" modules are plain text.

**Proposed clause** (to append inside A.2, in the item's own voice):

> A rendered Python-domain xref has TWO gates and a trap. `-W` catches nothing here; **`-n` catches
> every dead `:class:`/`:func:`/`:meth:`/`:attr:`/`:mod:` on an `.rst` page** (`[M]` 2026-09-21: 5
> of 5 under `-n`, 0 of 5 at default severity, two-sided control). What `-n` cannot see is the
> surface Sphinx never RENDERS: a docstring in an un-`automodule`'d module, and everything under
> `tests/`. The trap is `:noindex:`, which renders the docstring and mints NO target, so a role into
> such a module is plain text on a fully nitpicky build — `[M]` **24 of 48** `automodule` directives
> in this corpus carry it. check: `grep -c '^\.\. automodule::' docs/**/*.rst` for the population and
> read each directive's option block for `:noindex:`, then run `-n` as a pre-edit-vs-post-edit SET
> DIFF over the pages you touched — an absolute zero is unreachable while the plain-text convention
> stands. tell: a page whose roles all "resolve" by import and link to nothing.

**Founding case:** `[M]` 2026-09-21 probes, recorded in `archive_additions.md` §L-113, against
L-044 (the original, correctly-scoped `-n` measurement) and L-112 (the `:noindex:` two-sided HTML
control).

---

## U2. `vv-principles` — a published SENTENCE about which gates pin a claim is a coverage claim, with the same shelf life as a marker

**Joins:** `vv-principles` skill, § **"Log every caught bug"**, immediately after the paragraph
beginning *"A `catches(...)` or `verifies(...)` marker is a COVERAGE CLAIM with a shelf life, not a
topic tag"*; source `docs/development/skills/vv-principles.md`.

**Why.** The skill governs the MARKER and says nothing about the PROSE, yet prose is what future
readers quote: a theory page's "gates X and Y pin claim C", a close-out's residual table, an issue
comment's "this is covered by the sweep test". It decays by the same mechanism and is adjudicated
by the same evidence, but nothing currently says so — so the archivist, qa and test-architect each
re-derive it. The measured failure: a τ gate was credited in prose for reduced-operator arrays it
passes in 0.03 s under fully-garbaged factories, two screens after the same pass wrote the note
explaining that τ had LEFT that operator.

**Proposed clause:**

> The same holds for PROSE. A sentence claiming which gates pin a claim — in a theory page, a
> close-out, an issue comment or a report — is a coverage claim with the same shelf life as the
> marker, and it is adjudicated the same way: by a MUTATION that reddens the NAMED gates, never by
> topical adjacency. check: cite PER FIELD, not per topic (five arrays needed five different files;
> one had a sole catcher, another was cylindrical-only), and re-run the mutation whenever the gate's
> fixture, tolerance or budget moves. The highest-risk moment is REPLACING a gate you just demoted:
> the nearest-sounding sibling inherits neither its scope nor its teeth. tell: a prose citation of a
> gate whose docstring contains "until <phase> this asserted X"; a claim of coverage written in the
> same pass that explained why the quantity moved elsewhere.

**Founding case:** L-047 (the τ gate credited for arrays it cannot see) and L-102 (the ×2-pose gate
whose prose summary quoted the PRE-ruling assertion, twice, in two sections, while the gate's own
docstring recorded the date it inverted).

---

## U3. `instrument-doctrine` skill, X1 — a gate over a LIST must echo its input count, because an empty list and a clean tree print the same zero

**Joins:** the `instrument-doctrine` skill, § **X1**, under "Positive control, for a filter";
source `docs/development/skills/instrument-doctrine.md`.

**Why.** X1 already requires a known member of every shape before believing a filter's zero — that
covers a filter that is wrong. It does not cover a filter that is RIGHT running over an input list
that is EMPTY, which is a different mechanism with the same output and no member to validate
against. Measured twice in one session: an unsplit shell `$FILES` made two independent gates (an
xref probe and a markup scan) read ONE nonexistent path and print a clean `0`. The rule generalises
to every agent that writes a per-file gate loop, which is all of them.

**Proposed clause:**

> A gate that ranges over a LIST reports its INPUT COUNT beside its finding count, because a
> correct filter over an empty list and a correct filter over a clean tree print the same zero and
> the first one is a broken harness. check: print `len(inputs)` and assert it non-zero; assert one
> input resolves to an existing path. tell: a clean `0` from a loop whose inputs came from an
> unquoted shell variable, a `git diff --name-only` that matched nothing, or a glob with a typo.

**Founding case:** L-108 (both gates, one session, both flattering).

---

## U4. `code-search` — an UNQUOTED heredoc command-substitutes every backtick pair inside it, silently

**Joins:** the `code-search` rule, § **"`grep` here is ugrep, and an anchor inside an alternation
group matches nothing, silently"**, as a third mechanism beside the ugrep anchor and the
double-quoted-pattern case that rule already names; source `docs/development/rules/code-search.md`.

**Why.** `code-search` names zsh eating quotes and backticks out of a double-quoted pattern, and
notes that at least it prints `(eval): bad math expression` on a channel nobody reads.
`process-discipline` names the sibling for commit messages (`-F`, never `-m`). The heredoc case is
the one that is fully silent AND lands in the instrument itself rather than in a message: choosing
`<<PY` over `<<'PY'` so a path interpolates runs command substitution on every backtick pair in the
body, so a markup gate's own patterns lose their backticks before Python parses them. It failed loud
by luck (120 hits on clean prose); the same collapse in the other direction prints a clean 0, which
reads as a clean tree. Every agent writing a throwaway probe is exposed, and this project's probes
are full of backticks because its corpus is.

**Proposed clause:**

> A heredoc opened UNQUOTED (`<<PY`, chosen so a path interpolates) runs command substitution on
> every backtick pair in its body, so a pattern containing backticks is silently rewritten before
> the interpreter sees it. check: write every probe to a FILE under a QUOTED heredoc (`<<'PY'`) and
> pass paths in by environment variable; spell a backtick as `chr(96)`; and assert each pattern
> against its own witness before running it (`assert re.search(p, '<a known match>')`) — four
> one-line controls kill the whole class. tell: a gate that reports an implausible number on prose
> you just wrote, or a clean zero from a pattern you cannot see in the running process.

**Founding case:** L-102 (the gate whose four markup patterns collapsed to "match any bold"),
L-061, and its zsh sibling L-030 (an unquoted `$var` is not word-split, so a uniqueness loop ran
once on the concatenated string and printed a false 0).

---

## Considered and NOT proposed

- **"The brief is the FLOOR; live code is the rule"** — already `AGENT.md` Quality Checklist item 6
  for this agent, and `plan-authoring` §2 SHELF-LIFE / §4 RELAY for plans. A fifth copy is the
  duplication this pass exists to remove.
- **"Every universal carries its denominator"** — `instrument-doctrine` X2 and `plan-authoring` §2
  QUANTIFIER, with a dozen named shapes. The archivist's contribution is only that PROSE is one more
  population; the digest cites the clause instead.
- **"A refuted candidate is first-class output"** — `process-discipline` already carries it with the
  per-agent list.
- **A `-n` corpus-wide acceptance gate** — the obvious next step from U1, but it needs a corpus
  build to size its noise floor and this pass was forbidden one. Stated as the open question in
  §L-113 rather than proposed as a clause, because a rule clause prescribing an unmeasured procedure
  is the `plan-authoring` A-BRIEF'S-METHOD-IS-A-CLAIM defect.
