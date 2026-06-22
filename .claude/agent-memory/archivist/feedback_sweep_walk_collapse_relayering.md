---
name: sweep-walk-collapse re-layering docs
description: Updating SN sweep-architecture docs when a four-method direction×storage product (update_batch/residual_batch + apply_windowed/residual_windowed + SweepCellSlice packet) collapses to TWO storage walks × level-operation-objects × storage-free kernel pair — preserve Wave-2 history, add S6.4(e) current state, de-role retired symbols to literals, add a WHY-of-the-relayering note
type: feedback
---

When a within-module SN refactor COLLAPSES a direction×storage product surface
(the S6.4(e) graph-walk collapse, Issue #222) and you must update the theory
pages that documented the OLD surface, follow this shape. This is the sequel to
[[feedback_sn_internal_architectural_rewrite]] (the Wave-2 rewrite that
INTRODUCED the symbols this task retires).

**The collapse pattern in code terms** (read the code for ground truth — the
brief's prose is the map, the code is the territory): four direction×storage
walk methods (`apply`/`residual` full-field + `apply_windowed`/
`residual_windowed`) → TWO storage walks on the graph (`walk_full` oracle +
`walk_windowed` production), each parameterized by a LEVEL-OPERATION OBJECT
(`_CellSolve` solve-direction | `_CellResidual` apply-direction — direction is
an object, never a boolean). The per-level packet (`SweepCellSlice`) retired (it
only fed the retired storage adapters). The strategy's `update_batch`/
`residual_batch` → storage-free `cell_kernel_batch`/`residual_kernel_batch`
(PURE cell algebra, no gather/scatter — gather/scatter lifted to the walk layer,
ONCE, above every strategy).

**Rule**: update the "four primitives" table to the new factoring (typically:
OctantLabel / SweepDependencyGraph+its-two-walks / the level-operation pair /
the kernel pair), and ADD a load-bearing "WHY this re-layering" subsection with
its own `.. _<anchor>:` (here `sweep-dispatch-relayering`). The WHY is the
Cardinal-Rule-2 "build primitives not products" argument: a NEW closure would
have had to re-implement gather/scatter four times; lifting storage to the walk
layer lets a closure supply ONLY its kernel pair and inherit both storage
policies AND both directions for free. State the factoring explicitly: the
4-method product becomes 2(walks)×1(level-op pair, direction-by-object)×1(kernel
pair) where each factor varies independently.

**Why**: Cardinal Rule 3 — a future session reading the page alone must
understand both that the OLD four-method surface existed (history, so a `git
blame` / old-issue reader isn't confused) AND why the collapse is the principled
end state (so nobody re-introduces the product). The WHY-subsection is the
load-bearing content; the table is the index.

**How to apply** (the distinctive moves):
- **Preserve Wave-2 history via a `.. note::` "Architecture history — the
  dispatch surface re-layered twice."** Do NOT delete the Wave-2 narrative; flip
  it to past framing and add the S6.4(e) collapse as the current state. The
  retired symbol NAMES (`update_batch`/`residual_batch`/`SweepCellSlice`/
  `apply_windowed`/`residual_windowed`) stay as plain literals (double-backtick),
  explicitly tagged "appear below only as history".
- **De-role every retired symbol.** `:meth:`...update_batch`` → ``update_batch``
  (plain literal). An unresolvable `:meth:` renders as plain text WITHOUT a
  warning (Cardinal-Rule-1 staleness bug regardless), so `-W` will NOT catch a
  missed one — the only gate is the explicit grep
  `grep -n "update_batch\|residual_batch\|SweepCellSlice\|apply_windowed\|
  residual_windowed"`. After editing, re-run the grep and confirm EVERY
  remaining hit is an intentional historical literal (not an active role).
- **The Key Facts bullet** gets the new surface named (the two walks + the
  level-operation pair + the kernel pair) so the change is discoverable from the
  page header — same discipline as the Wave-2 note.
- **The performance follow-up note** (the Wave-2 "carry full-N buffers +
  octant_indices to cut levels×octants calls" deferral): preserve the historical
  direction but note that the SUBSEQUENT Phase 5 / S6.4 work took a DIFFERENT
  route to the same end (the rolling-frontier window). Don't pretend the old
  follow-up is still the open plan when a later wave superseded it.
- **Test-file rename**: `test_cell_update_batch.py` → `test_cell_kernel_batch.py`
  in the verification section; describe the new sha256 source-of-record pin on
  the two kernel bodies (the left-fold order is bit-identity-load-bearing). DROP
  hardcoded test counts (lesson L5) — say "term-level L0 on the kernel pair", not
  "10 tests".
- **The cross-doc changelog row** (index_convention.rst PR-INDEX-N table): a
  historical row naming a now-retired symbol (`DiamondDifference.update_batch`)
  is the symbol AS-OF-THAT-COMMIT. De-role to a literal and add a parenthetical
  forward-pointer to the rename: "``X`` (collapsed into :meth:`...cell_kernel_batch`
  at S6.4(e))". Preserves the changelog's historical accuracy while keeping the
  reader from chasing a dead symbol.
- **New active roles must resolve**: `walk_full`/`walk_windowed`/`cell_kernel_batch`/
  `residual_kernel_batch` are real methods (grep the code); `OctantLabel.streams`
  is a real property. `:ref:`wavefront-flux-cochain`` resolves cross-doc to
  operator_algebra.rst (grep `^\.\. _wavefront-flux-cochain:` across docs/ before
  citing).

**Build-gate note specific to this worktree**: a forced `-E` build surfaces a
WIDER baseline than the documented 1-warning (this run: 1 mesh.py paramref ERROR
+ 2 homogeneous.rst plotting WARNINGs + several verification.rst include-path
CRITICALs — the capability-matrix `.. include::` artifacts are an ENV artifact
needing a prior build, NOT staleness; see build-hygiene memory). NONE touch the
two edited files or the target symbols. The acceptance gate is "no
WARNING/ERROR/CRITICAL mentions discrete_ordinates / index_convention / any
target symbol", proven by grepping the build log for the filenames+symbols and
getting only progress-line matches (`reading sources... discrete_ordinates`),
never a severity line. The leftover `SyntaxWarning: invalid escape sequence`
lines are Python COLLECTING unrelated test files, not docutils.
