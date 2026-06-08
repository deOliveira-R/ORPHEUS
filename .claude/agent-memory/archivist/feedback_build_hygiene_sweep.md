---
name: build-hygiene-staleness-sweep
description: Pattern for a standalone Sphinx build-clean / staleness-fix task (not new theory) — forced -E rebuild to surface masked errors, module-move autodoc repoint, generated-artifact materialization, scope-guarded code-docstring flagging
metadata:
  type: feedback
---

Pattern for a docs-hygiene task whose goal is "make the Sphinx build
CLEAN + fix staleness" (NOT writing new theory pages), typically run
in a worktree with concurrent code work and a tight edit scope
(docs/ + conf.py only; flag any orpheus/ code-docstring fix).

(The forced-`-E` rebuild, the `WARNING:|ERROR:|CRITICAL:` grep, the
count-unchanged gate, and the venv/worktree facts are now in AGENT.md
"Build-Gating & Cross-Ref Reality" — this note carries only the
hygiene-task-specific fixes below. The `-E` discipline is ESPECIALLY
load-bearing here: one such task saw a plain build report 1 warning
when the real `-E` inventory was 13 errors + 9 warnings in unchanged
files.)

**Module-move autodoc repoint (the big one).** When `automodule::
orpheus.X.foo` fails to import, the module often MOVED *and* its
class surface changed. `orpheus.sn.quadrature` → re-exported
`orpheus.numerics.quadrature`; the five legacy subclasses
(AngularQuadrature / GaussLegendre1D / LebedevSphere /
LevelSymmetricSN / ProductQuadrature) collapsed into ONE `Quadrature`
value type with classmethod factories (gauss_legendre / lebedev /
level_symmetric / product). Verify the new surface with `python -c
"import …; print([dir])"` BEFORE writing refs — the class names are
dead, not just the module path. The doubly-stale `:class:` refs are
scattered across ~6 theory pages but are SILENT (nitpicky off → no
warning); only the `automodule` was build-breaking. Fix the
build-breaker; leave the silent theory-page refs unless explicitly
in scope (they render as plain code text, not broken links).

**DO NOT automodule a module with rich `:label:` docstrings.** Per
[[feedback_numerics_primitive_docs]]: pointing `automodule` at
`numerics.quadrature` surfaced "Inline strong start-string" +
"equation not found: <label>" warnings because submodules carry
`.. math:: :label:` blocks in docstrings (quadrature-selection-criterion,
product-mu-phi-cosines). The fix: replace the `automodule` block with
a `:class:`/`:meth:`/`:func:`/`:data:` cross-referenced prose section
(the same pattern numerics.rst uses for operator/measure/projection —
prose + refs, NOT automodule, "per-symbol docstrings live in the
modules, accessible via import path"). This eliminated both warnings
AND the refs resolve cleanly (0 warnings).

**Generated-artifact CRITICALs/plot-exceptions are ENV, not staleness.**
Two classes of "build problem" in a fresh worktree are missing
gitignored generated artifacts, NOT doc defects — DO NOT edit the
docs to "fix" them (that corrupts correct documentation):
  - `.. plot::` FileNotFoundError on `*.h5` → run the data converter
    `python -m orpheus.data.micro_xs.convert_gxs_to_hdf5` (reads
    tracked .GXS, writes gitignored .h5). ~10 large files, run in
    background.
  - `.. include:: ../_generated/*.rst` CRITICAL → run
    `python -m orpheus.derivations.generate_rst` (writes gitignored
    docs/_generated/). The main worktree HAS these; a fresh worktree
    doesn't.
  Confirm they're gitignored (`git check-ignore`) and that the
  generator only writes the artifact dir (no tracked-file change)
  BEFORE running. Materializing them is the principled clean-build
  fix; it leaves git status untouched.

**sphinx.ext.todo for method-implementer stubs.** `.. todo::` blocks
labelled "Archivist expansion needed" are the algebra-of-record
Sphinx-stub discipline — intentional rendered content, NOT to be
converted to `.. note::`. Enable `sphinx.ext.todo` + set
`todo_include_todos = True` so they stay visible to future sessions.

**Inconsistent-title-style fix = match the established level ladder.**
A `'` marker introduced as level-5 under an h3 parent whose OTHER
children use `^` (level 4) errors "skip from level N to N+2". Map the
file's marker ladder by first-appearance order (script: scan for
single-char underline rows), then convert the offending marker to the
sibling marker (`'`→`^`). Underline length must be ≥ title code-point
length (equal is OK). When the underline is short (e.g. 12 chars,
`'`-only), the old underline string isn't unique — edit with
title+underline as the match context.

**Malformed grid table = column-boundary drift.** A `:mod:` path that
was shortened left the cell's content shorter but the closing `|`
mis-padded → pipe lands at the wrong column. Script: list pipe/plus
indices per row; every row must share the border columns (e.g.
[2,10,34,80]); the offending row's last pipe is off. Recompute the
exact padding (interior width = next_border − this_border − 1) and
rewrite the one cell. Verify post-edit with the same index script.

**Scope-guarded code-docstring flag.** The standing baseline `:paramref:`
ERROR in `orpheus/geometry/mesh.py` (needs the uninstalled
`sphinx-paramlinks` ext) is a CODE docstring → out of a docs-only scope.
FLAG with the exact fix (`:paramref:`origin`` → `` `origin` ``, or
install+enable sphinx-paramlinks). Do NOT register a no-op role to
silence it — that produces a broken link, worse than the honest error.
The build is "clean apart from one flagged out-of-scope item".

**Stale `:doc:`/skills/X`` refs.** Skills are `.md` files outside the
Sphinx source tree — `:doc:` to them emits "unknown document". Convert
to plain-prose literal citations (`` the ``vv-principles`` skill
(``.claude/skills/vv-principles/``) ``) preserving the section-anchor
text. Matches the existing `` `.claude/scratch/…` `` literal style in
the same file.

**matrix.rst auto-gen drift.** `tools/verification/generate_matrix`
runs via a builder-inited hook on every build; the committed file had
4647, registry now 4760 (the ~113 drift). Run the generator (or just
build) to clear; it WILL re-stale when more tests land — the job is
clear the current drift + confirm pipeline health, not chase the
count.

**Quality scores (this task):** derivation-depth N/A (hygiene, not
theory); cross-references 5 (every quadrature factory linked);
code-traceability 5; the win was the FORCED -E rebuild discipline +
the env-vs-staleness triage that avoided "fixing" correct docs.
