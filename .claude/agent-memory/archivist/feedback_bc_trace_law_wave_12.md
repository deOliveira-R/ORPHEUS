---
name: BC trace-law Wave 12 docs synthesis (`docs/bc-trace-law-tensor-algebra-unification`)
description: Architectural-rewrite docs synthesis pattern for 12-wave refactor; 9-step wave-history → narrative compression; pyattr→attr/Greek-escape/eq-vs-ref label-vs-section disambiguation pitfalls; cross-method stub framing as scaffolding-not-mature-abstraction
type: feedback
---

# BC trace-law Wave 12 docs synthesis

Branch `docs/bc-trace-law-tensor-algebra-unification` 2026-05-11. Wave
12 of the 12-wave Boundary Operator refactor — the documentation wave
that consolidates Waves 0-11 narrative into durable Sphinx.

## Rule

A 12-wave refactor's documentation wave is the **load-bearing
deliverable**, not an afterthought. The code refactor IS the change;
the docs ARE the knowledge transfer. Treat it as MAXIMUM-EFFORT per
Cardinal Rule 3 — spend more time here than on any single code wave.

**Why:** Untracked closeout memos under
`.claude/agent-memory/method-implementer/` are per-session
artifacts. They are NOT durable. The Sphinx page IS durable. Pull
every salient narrative thread (design rationale, what-was-tried-and-
failed, ULP drift notes, semantic-correction motivation) INTO the
Sphinx page so the memos can be discarded without information loss.

**How to apply:** When summoned for a multi-wave refactor's
documentation wave:

1. Read the design plan in full (here:
   `.claude/plans/transient-giggling-cake.md` — ~1200 lines).
2. Read every per-wave closeout memo
   (`.claude/agent-memory/method-implementer/wave_*.md`).
3. Cross-reference with the Grand Report (here:
   `.claude/plans/neutron_transport_grand_report_v3.md` §16, §16A).
4. Build a master theory page (here:
   `docs/theory/boundary_conditions.rst`, ~1550 lines / ~7200
   words) that absorbs every salient piece.
5. Update the existing consumer pages (DO theory + operator_algebra)
   with a forward-pointer to the new master page + per-consumer
   detail tables.
6. Rewrite the package module docstring as the package's
   Sphinx-rendered architectural overview.

## Quality scores (this wave)

| Dimension              | Score | Notes                                                                                                                                   |
| ---------------------- | :---: | --------------------------------------------------------------------------------------------------------------------------------------- |
| Derivation depth       | 5     | Full §16A.3 three-layer decomposition + affine BC form + §16A.5 vacuum semantic correction.                                            |
| Cross-references       | 5     | 148 unique cross-refs in the new page; every concrete class, method, and module linked.                                                |
| Numerical evidence     | 5     | Wave-6 snapshot harness table with 8 cases × tolerance + ULP rationale.                                                                |
| Failed approaches      | 5     | 5-item anti-pattern catalog + the Option (a)/(b)/(c) Wave-7 retrospective for vacuum.                                                  |
| Code traceability      | 5     | Every realization branch in the SN dispatch table links to the realized primitive class.                                               |
| Derivation source      | 4     | No SymPy derivation script needed (this is architecture, not constructive math). The narrative inherits the Grand Report's structure. |

## Specific pitfalls learned this wave

### 1. `:pyattr:` is NOT a valid Sphinx role

Looks plausible (we have `:pyfunc:`, `:pyclass:`, …) but
`:pyattr:` is fake. Use `:attr:` for attribute references. Caught by
docutils ERROR "Unknown interpreted text role 'pyattr'" — these are
ERRORs not WARNINGs but `sphinx-build -W` does NOT promote them to
failure (it only catches Sphinx-level warnings, not docutils-level
errors). Run with `-E -W` and inspect "build finished with problems"
count to catch them.

### 2. Raw Greek `Γ_-` / `Γ_+` in text is parsed as RST reference

If `Γ_-` appears outside `:math:` / `\`\`...\`\`` quoting, docutils
parses it as a target reference and complains "Unknown target name:
'γ'". Fix: escape the underscore (`Γ\_-` / `Γ\_+`) or move the
expression inside a math role. The math role is preferred when the
context is mathematical; the escape is preferred for headers and
prose where math markup would be heavy.

### 3. `:label:` on a `.. math::` block is an EQUATION label

It is **NOT** a section anchor. Cross-reference equations with
`:eq:\`label\``, not `:ref:\`label\``. The previous-wave habit of
treating every `:label:` as a `:ref:` target produced one
"undefined label" warning that I fixed by switching to `:eq:`.
Conversely, section anchors are `.. _name:` directives placed
**above** the section title (with a blank line); they are reached
via `:ref:\`name\`` or `:doc:` / `:doc:\`<path>\``.

### 4. `:doc:` paths are RELATIVE TO THE SOURCE DIR

`:doc:\`../skills/vv-principles/error_catalog\`` is treated as a
docref to `<docs_source>/../skills/vv-principles/error_catalog.rst`
which doesn't exist. The fix is plain text:
"``.claude/skills/vv-principles/error_catalog.md``" wrapped in
literal-text role. Skill files are NOT in the Sphinx source tree.

### 5. `-W` exit code 1 from baseline ≠ regression

The baseline build had 9 pre-existing warnings. After my changes the
build still had 9 warnings — same set, same exit code. The
acceptance gate is "warning **count** unchanged pre/post-edit", NOT
"exit 0". This is the canonical close-out diff-style gate (per the
`feedback_retirement_docs.md` warning-count diff pattern); apply it
the same way for any docs change on a project whose baseline has
pre-existing warnings.

### 6. Cross-method stub framing matters

The user's main-agent augmentation explicitly steered the framing:
"don't over-promise the architecture". The 4 cross-method realizer
stubs (MoC / MC / CP / diffusion) are **scaffolding for future
modernization**, not a mature shared abstraction. The shared
`BoundaryRealizerBase` ABC is deferred per the "Unify after two
instances" rule until a second functional realizer ships. Document
the stubs as exactly that — the Protocol holds the dispatch table
in place; the abstraction emerges later. A `.. note::` admonition
block at the top of the stubs section makes the framing visible.

### 7. The vacuum semantic correction needs its own section

§16A.5 inflow-only-mask vs zeros-all is the single most subtle
design decision of the refactor. It needs its own section
(`bc-vacuum-semantic-correction`), the §16A.5 compatibility audit
preserved verbatim (13 production call sites enumerated), the
Wave-7 Option-(a)-over-(b)/(c) rationale captured, and inline
comments on the 3 retained-legacy test sites cross-linked.

## Forward-pointer issues filed

12 issues (#177-#188) covering all deferred follow-ups from the
12-wave brief. Each carries the wave-this-came-from reference + plan
file pointer so a future session can pick up the work without
re-reading the 12 closeout memos.

## What I would do differently next time

* **Build incrementally.** I wrote the full 1550-line master page in
  one shot, hit 9 docutils errors in the first build, and fixed them
  in a batch. Would have been cleaner to write the page in 3-4
  sections, build after each, fix errors locally before moving on.
* **Pre-check the LaTeX/RST role inventory.** Read `docs/conf.py`
  for the project's available macros (✓ I did this — `\Sigt`,
  `\Sigs`, `\keff`, `\kinf`, `\nSigf`) AND scan recent theory pages
  for the role conventions in use (I missed that the project uses
  `:attr:` not `:pyattr:`).
* **Diff the Sphinx warning count BEFORE the first edit.** Capture
  the baseline 9 warnings as the acceptance floor, NOT 0. Saves
  one round of confusion when the first `-W` build exits 1 for
  pre-existing reasons.

## Snapshot of deliverable sizes

| Deliverable                               | Lines      | Word count |
| ----------------------------------------- | ---------- | ---------- |
| New: `docs/theory/boundary_conditions.rst` | 1554       | 7210       |
| Update: `docs/theory/discrete_ordinates.rst` | +178/-48 | ~1300 new  |
| Update: `docs/theory/operator_algebra.rst`   | +165      | ~1100 new  |
| Update: `orpheus/geometry/boundary/__init__.py` | +302/-114 | ~2200 new |
| Update: `docs/theory/index.rst`               | +3/-1   | ~10 new    |

Sphinx build: 9 warnings (all pre-existing, unchanged from baseline),
exit 1 (baseline-equal). 148 cross-references in the master page.
6 equation labels with `:label:` tags (3 in new page, 3 in
operator_algebra additions). Nexus graph indexed the new page with
degree=144 (highly connected — every architectural concept
cross-linked).
