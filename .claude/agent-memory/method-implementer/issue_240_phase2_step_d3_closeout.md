---
name: issue-240-phase2-step-d3-closeout
description: #240 Phase 2 Step D3 — pure bit-identical mechanical rename CellUpdate*→DiscretizationScheme* and cell_update→scheme; macOS BSD-sed \b trap; worktree-safe git-grep scoping
metadata:
  type: project
---

# #240 Phase 2 Step D3 — the rename (BIT-IDENTICAL, NOT committed)

Branch `feature/sn-space-angle-tier2`. Host env → `.venv/bin/python`, `python -O -m pytest`.
Plan: `.claude/plans/issue_240_phase2_step_d_homing.md` §D3. NOT committed (user commits
after elegance + qa review).

**Why:** D3 resolves the `cell_update`-module-vs-attribute ambiguity and homes the renamed
generic advection–reaction discretization concept onto its native name. Pure rename, no
numerics. Followed D1 (`outgoing_face_from_average` + collapse closures) and D2 (home recon
ops onto the Base). NEXT = D4 (Σ-stateless realization), D5 (#239 2-D ScanMarch), D6 (docs).

**How to apply (future renames in this repo):** the two load-bearing traps below recur for
any large word-boundary rename on this macOS host with a sibling git worktree.

## The rename — Base-FIRST, word-boundary, perl (NOT BSD sed)

Substring collision `CellUpdate ⊂ CellUpdateBase ⊂ cell_update` forces Base-first ordering:
1. `CellUpdateBase` → `DiscretizationSchemeBase`
2. `CellUpdate` (Protocol) → `DiscretizationScheme` (Base already gone; `\b` spares the
   `IdentityCellUpdate`/`BadCellUpdate` test mocks)
3. `git mv orpheus/sn/spatial/cell_update.py orpheus/sn/spatial/scheme.py` (no collision)
4. `cell_update` → `scheme` (attr `mesh.scheme`/`self.scheme`, ctor param `scheme=`,
   `_CellSolve`/`_CellResidual` dataclass fields in `sweep_graph.py`, import paths,
   `:mod:`/`:class:`/`:meth:`/`:attr:` doc refs)

Plus consistency renames (no word-boundary so `\b` skipped them):
- test mocks `IdentityCellUpdate`/`BadCellUpdate` → `Identity/BadDiscretizationScheme`
- test names `TestCellUpdateBaseRegistry`→`TestDiscretizationSchemeBaseRegistry`,
  `TestDefaultCellUpdate`→`TestDefaultDiscretizationScheme`,
  `test_explicit_cell_update_honored`→`test_explicit_scheme_honored`,
  `test_cache_driven_sweep_matches_per_cell_update`→`..._per_cell_scheme_update`
- file `git mv tests/.../test_cell_update_protocol.py → test_discretization_scheme_protocol.py`
- doc path fix `discrete_ordinates.rst:1196` stale `tests/sn/spatial/test_cell_update_protocol.py`
  → correct `tests/sn/sweep/core/test_discretization_scheme_protocol.py` (was wrong-dir AND stale)

## ⚠️ TRAP 1 — macOS BSD sed does NOT support `\b`
The plan/brief recipe `sed -i '' 's/\b<old>\b/.../g'` SILENTLY does nothing on macOS (BSD sed).
Verified: a `\b` sed test left the line unchanged. **Use `perl -pi -e 's/\b...\b/.../g'`**
(PCRE `\b` works) OR BSD `sed -E 's/[[:<:]]...[[:>:]]/.../g'`. I used perl throughout.

## ⚠️ TRAP 2 — worktree safety via git-grep scoping
Sibling worktree `.claude/worktrees/nexus-workspace-wiring/` exists. `mcp__nexus__rename(apply)`
WOULD clobber it (it lists worktree files). Safety = `git grep -l` operates on the MAIN
checkout's tracked files ONLY and does NOT descend into a worktree (separate checkout) —
verified `git grep -l 'CellUpdate' -- '.claude/worktrees/'` returns nothing. NEVER filesystem
`grep -r`/`find`. Confirmed post-rename: zero worktree files in `git status`.

## Scope decision
Brief scoped to `orpheus tests docs`. EXTENDED to the tracked diag
`derivations/diagnostics/diag_phase_g_step2_cyl_si_fix.py` (it `from orpheus.sn.spatial.cell_update
import` — would break on the module mv; correctness > brief scope; NOT in forbidden set).
EXCLUDED `.claude/agent-memory/`, `.claude/plans/`, `.claude/lessons.md`, `.claude/scratch/`
(historical record — renaming history notes is wrong; brief's "history notes" = survivors).

## Registry keys SPARED (verified)
`"step"`, `"diamond_difference"`, `"linear_discontinuous"`, `"_no_kernel_strategy_test"` —
none contain the renamed tokens; `class DiamondDifference(DiscretizationSchemeBase,
key="diamond_difference")` keeps the string while the base renames.

## Intended residual SURVIVORS (git grep -n 'CellUpdate|cell_update' -- orpheus tests docs)
1. `discrete_ordinates.rst:2187` `test_cell_update_batch.py` — HISTORY NOTE ("S6.4(e) successor
   of test_cell_update_batch.py"); the file was migrated to test_cell_kernel_batch.py long ago.
2. `matrix.rst` — AUTO-GENERATED; regenerated via `python -m tools.verification.generate_matrix`
   → now reads `core/test_discretization_scheme_protocol` (5249 tests, count preserved).
3. `tests/_mutation/{README.md,diamond_spike.toml}` `test_cell_update_batch` — pre-existing
   mutation config referencing a long-gone file + a `field-role-typing` worktree path; not a
   `CellUpdate`-symbol token; out of D3 scope.
4. `test_cell_kernel_batch.py:22` — HISTORY NOTE.

## SHA gate — NO update needed
`test_cell_kernel_batch.py::TestKernelSourceOfRecord` 2 passed. The kernel method BODIES don't
reference `cell_update`/`CellUpdate` by name (only class headers / imports / docstrings do), so
`inspect.getsource` SHAs did NOT shift. (Anticipated by the brief.)

## Gate results — BIT-IDENTICAL ✅
- Clean import: `orpheus.sn.{solver,loss_representation,spatial.scheme,spatial.diamond}` OK.
- Residual `git grep -n 'CellUpdate|cell_update' -- orpheus tests docs` → only the 4 survivors above.
- Strict DriftWarning gate `tests/sn/sweep/core tests/sn/solve -W error::...DriftWarning`:
  **505 passed / 1 skipped / 4 xfailed** — IDENTICAL to pre-rename baseline (captured before
  touching anything). NO DriftWarning fired (the brief-anticipated vacuum_bulk_SLB 1-ULP warning
  did NOT appear at this HEAD — clean both pre and post).
- Route-around `tests/sn/{operators,spatial,sweep/core,sweep/cartesian_2d,solve}` -k route-around:
  **1083 passed** / 6 skipped / 7 deselected / 5 xfailed.
- `python -m tests._harness.audit` exit 0.
- Sphinx `-W` build exit 0 — all renamed `:mod:`/`:class:`/`:meth:`/`:attr:` refs resolve.
  Pre-existing `ld-cartesian-1d`/`ld-slab` Nexus skip-notices predate D3 (verifies()-label, not
  equation-node; 3 occurrences at HEAD).
- Diffs balanced (+46/-46 across diamond/geometry/loss_representation) — pure token substitution.

## NOT committed (user commits after elegance + qa review).
