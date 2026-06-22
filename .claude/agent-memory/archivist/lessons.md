# Archivist — Lessons

Read at the START of every invocation. Behavioral corrections only:
"what documentation/Sphinx/knowledge-architecture mistake did I make,
and what discipline did it teach?" The mechanical HOW lives in
`AGENT.md` ("Build-Gating & Cross-Ref Reality", "Close-Out Narrative
Arc") and in the preloaded skills (`vv-principles`, `algebra-of-record`,
`nexus-verification`). This file holds the PROCESS lessons those don't.
Campaign play-by-play (commit hashes, codenames) is retired — a lesson
here stands WITHOUT knowing what the campaign was.

The cross-cutting standard behind most of these: **a doc page is not
done when the prose reads well — it is done when every cross-ref
resolves against the LIVE tree, every claim's V&V level matches the
skill verbatim, every retired symbol leaves no dangling reference, and
the build's WARNING/ERROR/CRITICAL set is unchanged from the pre-edit
`-E` baseline.** Each lesson is one face of that.

---

## L-001 — Verify every claim against the LIVE code, NOT the quoted/docstring prose

The single most recurring trap. A task brief quotes "current stale
text", or a docstring describes a return shape, or a verdict memo
states a recommendation — and ALL THREE can be wrong relative to what
shipped. Three faces:

- **Brief quotes "stale text" that was already fixed.** The user's
  snapshot may pre-date a landed fix. Read the file (or `git show
  <prior-commit>:<path>`) before editing; if the current text already
  matches the desired "after", report "already fixed" rather than
  re-introducing the stale wording.
- **A docstring lies about a convention/shape.** A physics convention
  in a docstring (e.g. an index `r_i` vs `r_j`), or a return layout
  `(N, nx, 1, 1)` claimed while the array is `(N, ng, nx, ny)`, is
  ground-truth ONLY in the live code body. Read the array build / the
  consuming code; never transcribe the docstring's claim into a theory
  page.
- **A verdict memo records the RECOMMENDATION, not the OUTCOME.** An
  elegance/qa verdict's "BLOCKING NIT: drop X" may have been overruled
  — the shipped code took the alternative ("keep + strengthen test").
  Read the CODE to learn which resolution won.

How to apply: before citing any convention, shape, or design decision,
confirm it against the live source this session. Cross-refs that render
plain-text carry NO warning (see L-002), so a dead/stale `:func:` is a
Cardinal-Rule-1 correctness bug `-W` will never catch — grep the symbol.

---

## L-002 — Unresolvable code-xrefs render as PLAIN TEXT with no warning; this repo is NOT nitpicky

`-W` does NOT catch a dead `:func:`/`:class:`/`:meth:`/`:attr:` or a
stale alias-xref — they silently render as plain text. The acceptance
gate (count-unchanged from the `-E` baseline) only proves you added no
NEW warning; it is BLIND to staleness.

- After any carve that DELETES or RENAMES a symbol, `grep -rn "<symbol>"
  docs/` and repoint every hit on correctness grounds.
- Distinguish what DOES warn: undefined `[Key]_` citations warn (even
  non-nitpicky); intra-doc dangling `:ref:` warns under `-W`; cross-doc
  dangling `:ref:` renders plain-text. A new `:ref:` to a not-yet-
  existing section MUST create the labelled section in the SAME edit.
- Packages that are not member-`automodule`'d (transport.*, the
  operator/scheme/numerics leaves on several pages) render their
  `:class:`/`:meth:` as plain text BY PAGE CONVENTION. Match the page —
  do NOT half-surface 1–2 leaves by adding an `automodule` while the
  rest of the package stays plain. EXCEPTION: a new leaf whose
  docstring has NO `.. math:: :label:` is safe to automodule (no
  duplicate-label collision); a docstring WITH `:label:` must be
  cross-referenced in prose instead (automodule re-registers the label
  → duplicate-label / "equation not found" warnings).

How to apply: treat the grep gate (symbol exists / page-convention
matched) as the real cross-ref check; the warning count proves only
"added nothing new". (The forced-`-E` build, the three-severity grep,
and the venv/worktree facts live in AGENT.md — do not re-derive them.)

---

## L-003 — GREP the V&V matrix for an eq-label BEFORE renaming or removing it

An equation `:label:` may be a `@pytest.mark.verifies(...)` TARGET. If
it is, the test's oracle points at that exact name — renaming or
deleting it breaks the verification edge (and bumps a "no matching
equation node" line). The recurring mistake is rewriting a stale
equation's BODY and changing its LABEL in one motion, orphaning a
verifies-target.

- For a stale equation that IS a verifies-target: KEEP the label name,
  rewrite only the BODY (the claim is unchanged). Split a busy
  derivation into a label-preserving primary + NEW sub-labels for the
  decomposed steps.
- A label-rename ripples: an in-page symbol rename (`s_a → c_a`) flows
  to every `:eq:` and prose site that referenced it — a whole-page
  sweep, not a one-line edit.
- Section-label vs equation-label are different namespaces. `.. math::
  :label: X` → `id="equation-X"`, resolved by `:eq:`. A `.. _X:` anchor
  → `id="X"`, resolved by `:ref:`. When you need a section anchor and
  an equation sharing a name, suffix the section `X-section`. A
  `:label:` inside a CODE docstring is rendered by autodoc but is NOT a
  global `:eq:` target — cite the owning `:class:` and inline the math.

How to apply: grep the generated matrix (`docs/.../matrix.rst` or the
audit output) for every label you intend to touch FIRST; preserve
verifies-targets by name. Also `grep ':label:'` repo-wide before
MINTING a new label — duplicate labels across files collide (a real
warning), and a new label can collide with a same-named partition
predicate already living in another page.

---

## L-004 — Representational/structural eq-labels get a `.. vv-status: <label> documented` DIRECTIVE, not prose

A NEW equation that is a field-typing identity, a governing iteration,
a literature-transcribed definition, or a derivation-decomposition step
is NOT a solver claim. It must be tagged so the V&V matrix files it
under "Documented-only", not flagged as an unverified solver claim. The
recurring mistakes: (a) writing the status as prose instead of the
machine-read DIRECTIVE (a `--strict` audit then regresses); (b) leaving
a label that a NEW test's `verifies(...)` points at as an untagged
orphan.

- Structural/representational label → `.. vv-status: <label>
  documented` with a one-line rationale comment naming the bit-identity
  / foundation gate that pins the verifiable content.
- A label a NEW test `@pytest.mark.verifies(...)` points at is
  `implemented` (code + test, no eigenvalue/flux claim). If that test
  verifies a label that does not exist yet, CREATE the label in the new
  section — never leave a verifies-target orphaned.
- Pure derivation-decomposition labels (the affine-cell sub-steps, the
  facewise-separable tensor identity) sit untagged in a verified
  narrative chain by established page convention — match the page's
  siblings rather than inventing a status.

How to apply: for every eq-label you add, classify it (solver-claim /
representational / verifies-target) and apply the matching status
discipline. This is the V&V-vocabulary-curation duty (Directive 5) at
the label level — the matrix is the audit's source of truth.

---

## L-005 — A stub→rich expansion reads its sources in a fixed priority order; the docstrings are the prose SEED

Expanding a `.. todo:: Archivist` stub (the algebra-of-record handoff)
into rich narrative has a load-bearing source-reading order, and the
recurring mistake is writing from the brief alone:

1. **The close-out / verdict memo** — carries the bug ("confirmed live
   pre-fix"), the architecture-settled framing, the named tests, the
   verification numbers, the honest-interim state.
2. **The production docstrings** (the scheme/operator/class bodies, the
   `supports()` predicate, the SymPy `derive_*` docstrings) — the
   VERBATIM prose seed: the numerical-PDE statement, the contrast notes,
   the per-case table, the lit cites. These are the algebra-of-record.
3. **The test files** — the verification subsection: the named gates,
   what each asserts, the bit-identity-vs-principled-equivalence split.
4. **The SymPy module** (when present) — NEVER expand a stub without
   reading it; the narrative narrates it, does not compete with it. If
   you find an algebra error, return a DISPATCH_REQUEST for the
   method-implementer — never edit the SymPy yourself.

How to apply: read memo → docstrings → tests → SymPy, in that order,
before drafting. The honest scope (what shipped vs what is OWED to a
follow-on slice) comes from the memo's interim-state note — preserve it
verbatim; do not over-claim a wired-but-not-yet-iterating capability.

---

## L-006 — Cross-document duplicate citations are a real warning class; resolve, do not redefine

A `.. [Key]` bib entry duplicated across two standalone theory pages
emits a duplicate-citation warning. The recurring mistake is adding a
fresh `.. [Key]` block on a new page for a reference already defined
elsewhere.

- Before citing, `grep '^\.\. \[Key\]'` — if the entry exists on
  another page, cite it cross-doc (resolve), do NOT redefine.
- Where a page's existing convention already dodges the collision by
  citing a reference as PLAIN TEXT in the Literature list-table
  (because the `.. [Key]_` form would cross-doc-collide), MATCH that
  convention — add new literature as plain text too, no new bib entry.
- Pre-existing duplicate-citation warnings (cross-document cite
  collisions) are a known trade-off for standalone pages and do NOT
  need elimination during a close-out — verify the COUNT is unchanged
  pre/post, not that they are gone.
- FLAG (don't silently use) a conflated bib key when you spot one (e.g.
  a key whose title is one paper but whose cited equation content is
  another) — that is a Cardinal-Rule-1 correctness defect for the
  method-implementer/literature-researcher, surfaced from your
  cross-page vantage.

---

## L-007 — A retirement/relocation doc preserves the WHY and tombstones the wrong claim; it never deletes evidence

When an issue closes by FALSIFICATION (the approach cannot work) or a
type/decomposition is retired but its CONCEPT survives, the close-out's
archival value is HIGHER than a success story — it stops future
sessions re-attempting a dead path. The recurring mistakes are
rewriting history and deleting numerical evidence.

- **Preserve the motivation that LED to the investigation** — flip
  tenses ("is expected to" → "was expected to") but keep the logic. A
  future reader asking "why did anyone try this?" must find the answer.
- **Tombstone, don't delete.** When a new finding invalidates a
  published table/claim on the same page, add a `.. note:: **Retraction
  (date, Issue #N).**` immediately above it: (a) what the claim was, (b)
  why it's wrong (one line), (c) forward-pointer to the new section.
  Numerical VALUES stay; the INTERPRETATION gets the tombstone.
- **Retitle to the concept, KEEP the anchor.** When a type is retired
  but the concept survives in new realizations, retitle the section to
  the concept (not the dead type), KEEP the section anchor (cross-doc
  `:ref:`s auto-pick up the new title), and add a prominent succession
  note naming the new realizations. De-role dead module-path
  `:class:`/`:meth:` to literals (grep gate per L-002, not `-W`),
  past-tense the type, repoint present-tense claims to the realizations.
- The full 9-step CLOSED post-mortem arc + the PARTIAL/OPEN variant +
  the multi-issue audit-table pattern live in AGENT.md "Close-Out
  Narrative Arc" — follow it; don't re-derive it here.

---

## L-008 — Generated artifacts are NEVER hand-edited; materialize them, report the REAL number

Several Sphinx inputs are generated on every build by a `builder-inited`
hook (the V&V matrix via `generate_matrix`, the capability matrices via
`generate_capability_matrices`, the `docs/_generated/*.rst` includes via
`generate_rst`). The recurring mistakes are hand-editing the rendered
output and estimating a count instead of running the generator.

- A matrix/capability-table drift clears by running the generator (or
  just building) — never hand-edit the `.inc.rst`. Edit the registry-
  side metadata (`capability_rows()`) or the test that drives it.
- When a refactor changes how many rows the matrix auto-regenerates,
  REPORT the real post-regen number (it can differ from a brief's
  estimate — e.g. an auto-regen dropping 67→54 rows, not the est. −5).
  An auto-regen row delta is NOT a warning.
- In a FRESH worktree, missing generated artifacts (`.. plot::`
  FileNotFoundError on `*.h5`; `.. include::` CRITICAL on
  `_generated/`) are ENV gaps, NOT doc defects — materialize them
  (run the converter / generator, both write only gitignored dirs,
  confirm `git check-ignore` first). Do NOT "fix" the docs to route
  around a missing artifact — that corrupts correct documentation.
- Hard-coded test counts in RST go stale; `pytest --collect-only -q`
  for the current number rather than transcribing one.

---

## L-009 — Section-marker ladder is FILE-LOCAL; underline length is in CODE POINTS; reuse the file's existing same-depth marker

Recurring build-breakers from title machinery:

- The `=`/`-`/`~`/`^`/`'`/`"` underline ladder is per-FILE. Scan the
  file's first-appearance markers (grep the single-char underline rows,
  tally by depth) before picking a level, or you get "Inconsistent
  title style: skip from level N to N+2" — often at sections you did
  NOT touch (introducing a marker at a SHALLOWER depth than the file's
  existing one for that level pushes the old sections down a level).
- Underline length is measured in CODE POINTS, not bytes. An em-dash
  `—` or `÷`/`χ` is 1 code point but multiple bytes — size the
  underline with `len(title)` in Python, never `wc -c`. A global
  normalize-pass that touches PRE-EXISTING tolerated over-runs must
  RESTORE them (scope edits to your own lines).
- When adding a section at an existing depth, REUSE the file's existing
  same-depth marker character (grep for it). Introducing a different
  marker for the same logical depth is the classic "skip level" trigger.

How to apply: map the file's marker ladder first; size underlines with
`len()`; reuse existing depth markers. (The catalog of which file uses
which ladder is in AGENT.md — this lesson is the GENERAL discipline.)

---

## L-010 — V&V vocabulary in prose MUST match the skill verbatim; you are the curator (Directive 5)

You write the prose future readers quote when reasoning about
verification status. The recurring failure is paraphrasing a level
definition or over-claiming what a reference proves. The hard rules
(from `vv-principles`):

- **MMS does NOT verify eigenvalues** — it is source-driven, reaches
  flux-shape / convergence-order only. NEVER write "MMS verifies the
  eigenvalue". Eigenvalue claims need closed-form (`k_∞ = νΣ_f/Σ_a`) or
  semi-analytical references; SI≡Krylov twin agreement is
  necessary-but-NOT-sufficient (needs a structurally-independent leg).
- **L4 (code-to-code) proves zero correctness** — every L4 claim names
  its L0–L2 backing. NEVER "L4 proves correctness".
- **1-group eigenvalue is degenerate** (`k = νΣ_f/Σ_a` flux-shape
  independent) — NEVER "the 1-group test verifies the solver". (1G IS
  fine for a rate/convergence-order claim — declare the claim layer.)
- **Name the pillar** (closed-form / MMS / semi-analytical), not vaguely
  "analytical"; respect each pillar's evidence boundary.
- A **Mode-10 sub-floor term** is closed by STRUCTURAL teeth
  (producer-threads-at-machine-precision + consumed-sign-flip-moves-flux
  + a no-op control leg), NOT a tightened value band — and when no
  isolating regime exists (a boundary-trace slope below the bulk floor
  everywhere), say so: there is NO value-improvement leg to add. A
  prophylactic `.. warning::` ("do NOT write 'S9 recovers 2nd order at
  the boundary'") in the doc itself pre-empts the future over-claim.
- Bug attribution cites the failure mode (1–11) and matches
  `error_catalog.md`; a `catches(ERR-NNN)` claim is a coverage edge,
  not a topic tag.

**Skill-uplift duty:** when you notice a recurring published-prose
anti-pattern, a new failure-mode signature in a close-out, or a new
pillar-evidence-boundary case that `vv-principles`/`error_catalog.md`/
`algebra-of-record` does NOT yet capture, PROPOSE the skill edit in your
return. You read across all close-outs — the skill grows when you feed
it back. (Do not duplicate skill content into this file; point to it.)

---

## Quality self-assessment rubric (Directive 3)

Rate each output 1–5 and log the weakest dimension in the return:
Derivation depth · Cross-references · Numerical evidence · Failed
approaches · Code traceability · Derivation source (from `derivations/`,
never hand-written). The recurring weak dimension on TERMINOLOGY/ROUTING
fixes is "numerical evidence" — structurally absent (no flux moves → no
convergence table), not a deficit; say so rather than manufacturing one.
