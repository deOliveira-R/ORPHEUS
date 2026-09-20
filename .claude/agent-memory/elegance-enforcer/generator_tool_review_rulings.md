---
name: generator-tool-review-rulings
description: Reviewing a code/doc GENERATOR (docs -> .claude harness view) — the five probes that found every defect, and the "instrument with no reader" shape.
metadata:
  type: project
---

Reviewing `tools/docs/generate_harness.py` (W1-P3, 2026-09-20, branch
`docs/development-substrate`): a ~200-line tool copying `docs/development/`
pages into `.claude/`. Five probes produced every finding; reuse them on any
generator, linter, or `--check` tool.

**P1 — Does the check have a READER?** `ls -d .github`, `grep -rln <tool>
tests/`, and read how the build wires it. Here: no CI at all, no test, and the
Sphinx hook ran WRITE mode with `capture_output=True` discarding stdout — so
drift was silently *repaired*, never reported, and two docs asserted "`--check`
in CI". A `--check` whose unique arm (drift) has zero automated callers is
decoration. **This is the top finding shape for any generator.**

**P2 — Does the tool re-derive a quantity a LIBRARY in the same venv owns?**
The generator hand-wrote MyST's slug rule and extracted headings with a regex.
Diff its answer against the real parser over the whole corpus: 409 vs 396
anchors, 13 phantoms from `#` comments inside fenced code blocks; and the hand
slug drops non-ASCII word chars MyST keeps (`τ` -> `the--factor` vs
`the-τ-factor`). Note the grading: 0 live dead links, but the INSTRUMENT'S OWN
OUTPUT already differs — that is an X1 defect (teeth), not a corpus defect, and
it earns VIOLATION even though nothing is broken yet.

**P3 — Is the architecture's DIRECTION asserted by any line of code?** "The
flow is one-way, docs -> .claude" was a manifest naming convention: nothing
constrained `target` to `.claude/`, so a copy-paste typo would overwrite the
generator's own source in write mode. Ask of every stated invariant: which line
fails if it is false?

**P4 — Instrument branch/arm activation counts.** Instrument the run, don't
read it. Measured: `relative_to` try-arm 0 of 246 (dead; the fallback
re-implements `os.path.relpath`, agreeing 246 of 246); budget `None` arm 0 of 9
while the docstring said "EVERY text target has a budget" (anti-#20); 236
`headings()` calls over 26 files.

**P5 — Orphan/coverage census in BOTH directions.** Manifest sources vs docs
pages; stamped outputs vs manifest targets. Both were 0 here, but note the
asymmetry: a source missing from the manifest is caught by Sphinx's
orphan-toctree warning, a REMOVED manifest entry leaves a stamped `.claude/`
file loading into every agent forever with no instrument at all.

**Grading ruling — do NOT flag `kind` strings as stringly-typed dispatch when
they are discriminated ONCE.** rule/skill/index/agent were enumerated at one
site and branched at one site. That is the skill's positive form
("discriminate once, at the boundary"), not anti-pattern #4. The real (NIT)
finding was that two of the four names share one code path and the contract
page documents only three of them.

**Derived-constant tell worth carrying:** `begin_line = BEGIN.split(" — ")[0]`
— a stable marker prefix obtained by parsing the tool's own format string. Edit
the em-dash and every marker lookup misses, so write mode APPENDS a duplicate
generated block to nine files instead of replacing. Name the prefix; compose
the template from it.
