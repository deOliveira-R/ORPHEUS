# `AGENT.md` — corrections and promotion candidates

Source of truth: `docs/development/agents/archivist.md` generates only the role block at the top;
the body of `.claude/agents/archivist/AGENT.md` is hand-maintained, so these are edits to that file
unless the harness says otherwise. I edit nothing — this is the proposal.

Two **CORRECTIONS** first, because a false clause in the always-loaded definition is a
Cardinal-Rule-1 bug and it has been deciding real rename verdicts. Then two **PROMOTIONS**, kept
deliberately small: a thin AGENT.md of sharp identity principles beats a bloated one, and the
promotions below sharpen clauses that already exist rather than adding sections.

---

## C1 (CORRECTION, blocking) — the cross-doc `:ref:` claim is FALSE, and the `-n` claim is over-broad

Two places assert it, and they are one source, not two — one copied the other
(`instrument-doctrine` X4).

**Site 1, Quality Checklist item 3**, last parenthesis:

> "(Undefined `[Key]_` citations and intra-doc dangling `:ref:` DO warn; cross-doc `:ref:` renders
> plain-text.)"

**Site 2, "Cross-ref reality (this project is NOT `-n` nitpicky)"**, bullets 1 and 3:

> "Unresolvable `:func:`/`:class:`/`:meth:` refs render as PLAIN TEXT with NO warning — **and `-n`
> (nitpicky) does NOT save you either.** … Intra-doc dangling `:ref:` IS caught by `-W`; cross-doc
> dangling renders plain-text."

**`[M]` 2026-09-21** (two throwaway Sphinx projects, this repo's venv, live targets beside dead
ones; full transcript in `archive_additions.md` §L-113):

- A dangling `:ref:` warns at **DEFAULT severity** whether the target is intra-doc or cross-doc, and
  whether the role is bare or carries explicit text — three `[ref.ref]` warnings, identical subtype.
  A dangling `:doc:` warns as `[ref.doc]`. **"cross-doc `:ref:` renders plain-text" is false.**
- Five dead python-domain roles on an `.rst` page: **0 of 5** warn at default severity, **5 of 5**
  warn under `-n`. So "`-n` does not save you" is true of the surface Sphinx never RENDERS (a
  docstring in an un-`automodule`'d module; everything under `tests/`) and **false of the `.rst`
  corpus**, which is rendered by construction — and the `.rst` corpus is where the archivist works.

**Why it matters, not just that it is wrong.** L-063 KEPT a stale eq-label prefix because "a
cross-doc dangling `:ref:` renders plain text at every severity — renaming buys cosmetics and risks
a silent break", over 8 cross-doc citers; L-076 renamed an anchor because "L-063's
silent-cross-doc-break caution did not bind" at 1 citer. A false premise was setting rename policy
in both directions. The KEEP verdicts survive on their real grounds (a stale NAME is not a false
CLAIM; many citers is a cost) but the silence argument must be struck.

**Proposed replacement, item 3's parenthesis:**

> (Undefined `[Key]_` citations **and every dangling `:ref:` / `:doc:`, cross-doc included, DO
> warn** at default severity — `[M]` 2026-09-21 — so both label classes are rename-gated by the
> build. What is silent is the code-xref.)

**Proposed replacement, "Cross-ref reality" bullet 1's flag clause** (keep the rest of the bullet
verbatim — the grep gate, the `docstring of <module>` probe and the Cardinal-Rule-1 framing are all
correct):

> **Unresolvable `:func:`/`:class:`/`:meth:`/`:attr:`/`:mod:` refs render as PLAIN TEXT with NO
> warning at default severity — but `-n` DOES catch them on an `.rst` page** (`[M]` 2026-09-21: 5 of
> 5 under `-n`, 0 of 5 by default). What no severity can see is the surface Sphinx never RENDERS: a
> docstring in an un-`automodule`'d module, and every file under `tests/`. **And `:noindex:` is a
> second plain-text mechanism** — it renders the docstring and mints no target — `[M]` 2026-09-21
> **24 of 48** `automodule` directives in the source carry it. Run `-n` as a
> pre-edit-vs-post-edit SET DIFF over the pages you touched; an absolute zero is unreachable while
> the plain-text convention stands. The `grep` gate is still the acceptance evidence, for the
> reason below.

**Delete bullet 3** ("Intra-doc dangling `:ref:` IS caught by `-W`; cross-doc dangling renders
plain-text") and keep only its imperative, which is right: *when you introduce a `:ref:` to a
not-yet-existing section, create the labelled section in the SAME edit — the build will otherwise
fail, not rot.*

Also propose retitling the section: **"Cross-ref reality (this project does not RUN `-n` — but `-n`
sees the `.rst` corpus)"**. The current title asserts the refuted generalisation in the one line a
skimming reader keeps.

---

## C2 (CORRECTION) — AGENT.md quotes a frozen count, in the section that forbids quoting frozen counts

Same section, bullet 1: *"the doc source carries only ~45 live `automodule` directives
(2026-08-03)"*. `[M]` 2026-09-21: **48**, of which **24** are `:noindex:` — so the number drifted
and, worse, it measures the wrong thing, because half of the "surfaced" modules mint no targets.

This is the section that two paragraphs earlier says *"a quoted baseline is a frozen claim that
rots exactly like a stale `NEXT =` pointer: it reads as authority long after it stopped being
true."* The fix is the one that section already prescribes for the warning baseline: **publish the
command, not the number.**

> **check:** the population is `grep -c '^\.\. automodule::' docs/**/*.rst` (excluding `_build`),
> and the plain-text half is however many of those carry `:noindex:` in their option block — read
> it, do not quote it. `[M]` 2026-09-21 it was 48 and 24; if you are citing that pair, re-measure
> first.

While there: the same section's list of packages *"NOT automodule'd anywhere"* (`transport.*`,
`numerics.spaces.*`, module-level private `_helpers`) is a roster, and a roster is a universal
owing its denominator (`instrument-doctrine` X2). It was not re-measured this session. Propose
adding *"re-derive this roster before relying on it; it is a snapshot"* rather than extending it.

---

## P1 (PROMOTION) — sharpen Quality Checklist item 6 with the two clauses that make it enforceable

Item 6 already carries the identity principle ("Verify EVERY claim against the LIVE source this
session — not the brief, the docstring, or a verdict memo") and it ends *"confirm cited numerical
results are reproducible."* That last clause is the one that keeps failing, because it is read as a
spot check. Two sentences make it a procedure, and both are applied on essentially every task:

> Re-derive the pass's numeric literals as a **SET, in one script, at the END** — a spot check finds
> neither of the two defects the set finds (a figure that was invented, and a figure that was real
> and answered a different question). And when two honest measurements disagree, do not adjudicate:
> find the PREDICATE or the STATISTIC that makes both true and publish the arithmetic, because an
> exact reconciliation proves the census complete and names what it excluded.

Founding cases L-109 (both defects in one run) and L-107/L-108/L-100 (the reconciliations). The
digest keeps the instances; item 6 gains the standing bar.

---

## P2 (PROMOTION) — the two-build session sequence belongs in the definition, not in recalled lessons

There is no "how a docs session is sequenced" statement anywhere in AGENT.md: the Build-Gating
section covers what to run, never in what ORDER, and the consequence is measured — the same pass
cost four builds once and five another time, every extra one bought by an edit made after launching
a verification build. It is identity-level (it shapes every task's shape), it is one sentence, and
it has no home in any rule or skill.

Propose it as a new numbered item at the head of the Quality Checklist, before item 1, since it
orders the items that follow:

> 0. **Sequence the session so you build TWICE:** baseline `-E -W` → *all* edits → *all* residual
>    greps → the xref gate and your own import probe → the AST doc-only proof → ONE verification
>    build. The self-consistency pass over prose YOU authored — universals, quotations,
>    denominators, superlatives, symbol collisions, aspirational rows — runs to EXHAUSTION *before*
>    the first verification build, never interleaved with it. A residual grep always finds one more
>    site; every build launched before it does is wall-clock spent for nothing.

Founding cases L-054, L-064, L-081.

---

## Considered and NOT promoted

- **"THE SPINE — a page is DONE when …"** — it is the Quality Checklist's own bar in one line, and
  the checklist is already always-loaded. It stays at the head of the digest as the reading
  shortcut, pointing at the checklist; promoting it would put the same bar in two places, which is
  what this pass is undoing. (One clause of it IS worth landing: checklist item 2 says "Zero
  warnings", while the section 300 lines below correctly says the gate is the count-and-SET diff
  from a freshly-measured `-E` baseline. Item 2 should point at that section instead of stating a
  different bar.)
- **The nine-step close-out arc, the build-gating mechanics, the venv/worktree facts** — already
  AGENT.md's, correctly.
- **The xref-gate role-scoped blindness (digest §2b)** — a measured lesson with two controls and an
  unlanded one-line fix in `tools/check_docstring_xrefs.py`. It belongs in the digest while the fix
  is pending; it becomes an AGENT.md line only if the fix is abandoned. Flagging it as the one item
  whose home should be revisited after the gate is patched.
- **Every §6 event-class shape** — 30 recipes; they are exactly what a digest is for.
