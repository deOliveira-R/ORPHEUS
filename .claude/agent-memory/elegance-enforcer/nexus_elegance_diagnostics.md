---
name: nexus-elegance-diagnostics
description: Five Nexus structural-smell tools (native_place/twin_paths/discriminations/dead_functions/protocol_conformers) that mechanize my review axes 1-4+7; how to route them and their FP classes.
metadata:
  type: reference
---

A Nexus tool family (shipped in sphinxcontrib-nexus branch `feat/native-place-diagnostic`,
2026-06-17) mechanizes the structural-smell first pass of an elegance review. Each is a
`GraphQuery` method + MCP tool + `nexus <cmd>` CLI; all read-only, all take
`exclude=""` (comma-sep substrings on top of an is_test filter) + `limit`, all emit JSON
with `file_path`/`lineno` per node (filter to a diff by path). Verify the branch merged
before relying on these in a session: `nexus <cmd> --help`.

Drafted a NEW skill `nexus-elegance` for these (the PRE-refactor discovery phase; sibling to
`nexus-refactoring` which EXECUTES). Decision rationale: refactoring skill is rename-only
(WHAT to change); these find WHERE. Verification skill is Rule-1 correctness; these are
Rule-2 architecture. One home, not split — all five share "candidate-not-verdict + known FP
class". Main agent owns the file writes.

**Tool → my review axis / coding-elegance pattern:**
- `twin_paths(min_similarity,min_tokens)` → axis 2 path-multiplicity / Pattern 2 / anti-#1.
  Catches array-math clones the CALL GRAPH can't see (einsum, slicing). FP: symmetric-by-design
  (apply/transpose, solve/residual = principled mirroring, my `_CellSolve`/`_CellResidual`
  ruling); test fixtures + student_resources teaching copies fire at 1.0 → exclude mandatory.
- `native_place(min_callers)` → coding-standards "native place" / Feature-Envy / supports Pattern 5.
  Free fn whose every caller is ONE class → move in. FP: `likely_free_primitive=true` = drop
  (correctly free, my `_initial_guess_values` pure-extractor ruling).
- `discriminations(min_sites)` → anti-#3/#4 (boolean/stringly dispatch) via new `discriminates_on`
  AST edge. ≥2 sites branching one tag → polymorphism candidate. FP: single dispatching site
  IS the centralization (my "one shared predicate=STRENGTH, a SECOND spelling=smell" ruling).
- `dead_functions` → axis 7 dead weight / aggressive-retirement / anti-#11/#12. Zero non-test
  `calls` in-edges. CANDIDATE only — dynamic dispatch invisible. Trust private+undecorated;
  decorated=registry hook, public=likely API (both FP-prone).
- `protocol_conformers(min_methods)` → Pattern 4 (declare contract → checker enforces) +
  Pattern 6 (≥2 conformers earns the declaration). NAME heuristic; pyright authoritative.

**Workflow:** PR review = run whole-graph, filter output to changed-file set by `file_path`.
Debt sweep = run unfiltered → Architectural Opportunities + Rule-4 issues. Default
`exclude=student_resources,derivations,scratch,tools` (without it twin_paths/dead_functions
top results are student_resources teaching dupes). Every hit is a QUESTION; apply the
bug-habitat test by reading before promoting to VIOLATION.

**Live-graph corroboration signal worth keeping:** `native_place` AND `discriminations` both
independently flagged `axes_from_legacy_mesh` and `_subdivide_zone` on `/tmp/orpheus_full.db`.
Two tools agreeing on one symbol (homeless fn that ALSO branches a repeated tag) is a stronger
escalation than either alone. The `coord` discrimination was 16 sites live (the standing
geometry-tag cluster I already track).
