---
name: canonical-convention-page
description: Pattern for shipping a multi-PR layout/storage migration's docs deliverable as a single canonical Sphinx theory page (`index_convention.rst` for Issue #196 PR-INDEX-6 — SN principled storage flip). Distinguishes the "canonical statement" page from per-PR closeout memos and from the migration plan; covers the 9-section anatomy (Key Facts / Overview / Derivation / Cross-section / History / Step-1 gate / Numerical evidence / Layout reference table / Gotchas / Future work / Cross-references / References); sweep-audit classification rubric for keep-vs-flip decisions on legacy-shape mentions.
metadata:
  type: feedback
---

# Canonical-convention page (Issue #196 PR-INDEX-6 anchor)

A multi-PR architectural migration's docs deliverable is a single
**canonical Sphinx page** that consolidates the convention into one
quotable, future-readable artifact.  This rule is the load-bearing
output of the final-cleanup PR in a migration sequence — the one
that retires bridge-transposes from code and leaves the
documentation-of-record behind.

**Why:** The 9-step closeout-narrative arc (in MEMORY.md) is for
ISSUES.  The canonical-convention page is for **theory docs**.  A
closeout memo lives at `.claude/agent-memory/`; a canonical-
convention page lives at `docs/theory/`.  The two serve different
audiences: the closeout is what the orchestrator/main-agent reads
to gate-keep a PR; the canonical page is what every future session
quotes when reasoning about the convention.

**How to apply:** When dispatched to ship docs for the
**final-cleanup PR** of a multi-PR refactor (a "PR-N+1 docs PR"
after PRs 1..N flipped semantics), produce ONE rich theory page
with this 13-section anatomy.

## The 13-section canonical-convention page anatomy

1. **Title + `:contents:` directive.**  Title carries the
   convention itself in the headline (e.g., "SN Index Convention --- ``(N, ng, nx, ny)``").
2. **Key Facts admonition.**  Authoritative statement of the
   convention.  Bullet list of: the shape rule, the orthogonal
   shape rule (e.g., scalar vs angular), the singleton-axis rule
   for the degenerate-dimension case, the cell-flattening (or
   analogous) invariant, the four-operator algebra contract.
   Close with an "authoritative-origin" admonition pointing at the
   §Derivation + §History + §Numerical-evidence sections.
3. **Overview.**  Frame the four (or however many) indices,
   contrast with the legacy convention, name the migration sequence
   length.
4. **Derivation --- why this convention is principled.**  The
   mathematical reasoning.  Full prose elaboration of the
   axis-priority table from the migration plan §1.  Include a
   labelled `:eq:` equation if the convention follows from an
   operator-algebra structure (e.g., within-group block-diagonality).
   Add an "Algorithmic consequence" subsection showing the win.
   Cite the textbook source ([LewisMiller1984]_ etc.).
5. **Companion convention** (e.g., cross-section storage paired
   with flux storage).  Includes a labelled `:eq:` for the
   round-trip identity if there is one (e.g., cell-flattening
   `legacy[i,j,g] == new[g,i,j]`).  Code excerpt of the
   `__debug__` assertion that pins the identity.
6. **History --- the N-PR migration.**
   - **The proposal that was wrong.**  Quote the discarded design
     (often from an earlier memo); explain why it was wrong; cite
     the principle that flagged it (`coding-elegance` Pattern N).
   - **The N PRs.**  list-table with commit hash, scope per PR.
     This is the canonical place to record "PR-INDEX-3 lives at
     commit X" — readers can `git show` that commit.
   - **What stayed deliberately legacy** subsection.  Document any
     deferral (e.g., FD-matvec internal contract) with the scope
     and the forward-pointer to the resume PR.
7. **The load-bearing verification gate.**  The verbatim
   bit-identity check or analogue.  Include the actual Python from
   the most-load-bearing PR's closeout memo (paste verbatim — don't
   paraphrase).  Quote the max-residual number ("max abs diff
   1.75e-14") because the reader will Google it.
8. **Numerical evidence.**
   - Regression / verification snapshot list with residual
     list-table (verbatim numbers from the most-recent closeout).
   - L0 / L1 gate inventory (pre / post).
   - Performance benchmark (mean ms or whatever the unit is) as a
     list-table showing the per-PR progression.
   - Equivalence-class checks (`nulp=...` for principled-equivalence
     cases).
9. **Layout-by-array reference table.**  list-table mapping every
   array a future maintainer encounters to its shape + definition
   site.  Include the "exceptions" subsection (apparent layouts
   that aren't actually exceptions; e.g., primitive contracts like
   the `(nx, K, ng)` scan-axis-leading contract).
10. **Gotchas and subtleties.**  Singleton-axis rule, sibling-
    convention reminders (e.g., SigS[g_from, g_to] unchanged),
    per-material vs per-cell distinctions, test fixture order.
11. **Future work.**  Each deferral gets a labelled subsection
    with: PR-Nth scope estimate (~XXX LoC), bullet list of files
    touched, dispatch sub-agent (method-implementer for code,
    archivist for docs).  Cite the migration plan's §7 (or
    analogous) deferred-work register.
12. **Cross-references.**  Bullet list of every adjacent page,
    plan file, closeout memo.
13. **References.**  `[Citation]_` blocks with section numbers
    cited inline.

## Sweep-audit classification rubric

The companion deliverable to the canonical page is a sweep audit:
grep for legacy-convention mentions across `orpheus/` and `tests/`,
classify each, flip the ones that can be flipped.  Use this rubric:

| Class | Description | Decision rule |
|---|---|---|
| **Docstring describing the convention** | Plain shape note in a `"""..."""` block | **FLIP** to principled + cross-ref to canonical page |
| **Comment naming the legacy convention as part of an invariant check** | A `__debug__` block, a round-trip identity comment | **KEEP** — the comment NAMES the legacy convention because the identity tests legacy↔principled |
| **Primitive-contract description** | A scan-axis-leading description, a Blelloch leading-axis comment | **KEEP** — primitive contracts (per Pattern 7) are local to their definition site |
| **Documented intentional construction adapter** | A "build legacy, view-transpose principled" pattern preserved for snapshot bit-identity | **KEEP** — flipping the construction order changes numerical values |
| **Pre-refactor narrative in a hand-invoked diagnostic** | A structural-audit text block describing legacy code | **KEEP** — diagnostics are hand-invoked; partial flips corrupt the audit's meaning |
| **Generic primitive test using SN-leaning variable names** | A test of a generic primitive (e.g., `HarmonicMomentProjection`) that broadcasts across any trailing axes, but happens to call them `(nx, ny, ng)` | **RENAME** to abstract names (`a, b, c`) + add a docstring stating the primitive is layout-agnostic |

For the FLIP class, the rewrite preserves the shape note (don't
delete it — replace), and ADDS a `:ref:\`theory-canonical-page\``
cross-reference.

## Code docstring-flipping idiom

The "PR-N era docstring" pattern (a docstring written during PR-N
that describes bridge-behaviour) needs the present-tense flip the
final-cleanup PR provides:

```rst
OLD ("will retire"):
    PR-INDEX-1 carries the principled (N, ng, nx, ny) internally;
    public signature is unchanged — entry/exit transposes convert
    caller-side (nx, ny, ng) inputs to principled layout.

NEW ("post-migration"):
    PR-INDEX-1 through PR-INDEX-5: internal arrays AND the public
    signature both carry principled (N, ng, nx, ny) (see
    :ref:`theory-canonical-page`). No entry/exit transposes are
    required — caller-side principled inputs flow directly.
```

The flip is: (1) extend the PR-range, (2) collapse the future-tense
"will retire" to present-tense "are gone", (3) add the canonical-
page cross-reference.

## Sphinx-build acceptance gate

The acceptance gate is **baseline-warnings-unchanged** (per the
`feedback_bc_trace_law_wave_12` rule), NOT count=0.  Run
`sphinx-build` pre-edit + post-edit and verify the warning count is
identical AND that no NEW warning text appears.  The new
canonical-convention page must introduce ZERO new warnings.

If `:doc:` cross-references emit `unknown document`, the fix is
ALWAYS to use the actual file path (`docs/api/discrete_ordinates`,
NOT `docs/api/sn`).  Don't invent doc paths to match cleaner names.

## Page-length budget

Aim for **800-1000 lines** of RST.  The "Sphinx-as-brain" payoff
beats the "concise summary" trade.  Compression-driven cuts kill
the page's ability to serve as a standalone expert briefing.

For the issue-196 PR-INDEX-6 instance: 931 lines, 107 KB HTML
output, 13 sections, 7 list-tables, 3 labelled `:eq:` blocks, ~12
cross-references to other theory pages, 3 citations.
