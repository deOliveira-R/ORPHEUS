# Nexus bug report — `retest --scope=branch` returns 0 changed symbols

**Date**: 2026-05-14
**Repo**: `ORPHEUS`, branch `refactor/sn-operator-algebra`
**Tool**: `sphinxcontrib-nexus` MCP server
**Severity**: Blocks dev-workflow integration of `nexus.retest` for diff-aware test selection.

---

## Symptom

`mcp__nexus__retest --scope=branch` returns `{must_retest: [], should_retest: [], changed_symbols: [], total_tests: 2170, safe_to_skip: 2170}` despite the branch having **179 commits ahead of `main`** with **89,001 insertions / 2,107 deletions across 296 files**.

Likewise, `mcp__nexus__detect_changes --scope=branch` returns `{changed_symbols: [], affected_symbols: [], total_changed: 0, total_affected: 0}`.

Confirmed via:

```bash
git -C /Users/rodrigo/git/nuclear/ORPHEUS rev-list --count main..HEAD
# → 179

git -C /Users/rodrigo/git/nuclear/ORPHEUS diff --stat main..HEAD | tail -1
# → 296 files changed, 89001 insertions(+), 2107 deletions(-)
```

Modified files include core Python symbols (`orpheus/sn/sweep.py`, `orpheus/sn/spatial/sweep_cache.py` NEW, `orpheus/sn/spatial/diamond.py`, `orpheus/sn/geometry.py`, plus ~290 others) that should appear as `changed_symbols`.

## What works

- `--scope=unstaged` returns `0` correctly — working-tree modifications are only `.mcp.json` and `docs/verification/matrix.rst`, neither containing Python symbols.
- `--scope=staged` (untested but plausibly similar to unstaged in semantics).
- `--scope=all` returns the same `0` as `branch`, suggesting the same git-diff base resolution path.
- The underlying graph is healthy:
  - `nexus.stats` → 9600 nodes, 85013 edges, 0 stale docs / 114 checked.
  - `nexus.verification_audit` → 2442 explicit `verifies(...)` test-edges + 27324 call-graph-inferred test-edges, 444 equations (376 verified, 58 documented, 10 implemented-but-untested).
- The graph reflects post-2.6-Q3 state correctly (e.g., `iter_cell_visits` / `iter_cells_by_direction` no longer in the graph; `dag_walk` is the sole iteration method).

## Hypothesis

The `branch` and `all` scopes' git-diff base resolution is broken. Three candidates:

1. **The git-diff base is hard-coded to `HEAD`** instead of `main` (or the configured base branch). HEAD vs HEAD = 0.
2. **The branch base is read from a config that's stale or unset**, defaulting to HEAD.
3. **Symbol mapping happens post-rebuild, and the graph indexes HEAD only**, so even with a correct git diff there's no "previous" symbol set to compare against. The graph would need a snapshot of the pre-diff symbol IDs.

Hypothesis 3 is the most architecturally interesting — it'd mean `retest` needs either (a) a stored snapshot of `main`'s symbol IDs at last rebuild, or (b) two graph instances (main + HEAD) to diff.

## Reproduction

In ORPHEUS at branch tip `01999e0`:

```python
# MCP call: mcp__nexus__retest, args: {"scope": "branch"}
# Expected: must_retest contains tests for changed Python symbols (sweep_cache, diamond, sweep, geometry, operator, solver, ...).
# Actual: empty.
```

## Impact

Blocks **Tier B slow-test selection** for the ORPHEUS SN suite. The intended workflow:

```bash
TESTS=$(nexus retest --scope=branch | jq -r '.must_retest | join(" ")')
pytest $TESTS
```

would replace `@pytest.mark.slow` manual tiering with graph-derived selection. Currently impossible because retest returns 0 must-retest tests for any diff scope larger than unstaged.

## Verification graph (for reference)

The Nexus graph at branch tip `01999e0` has the right shape:

- `tests_declared: 2442` — explicit `@pytest.mark.verifies(...)` decorator markers.
- `tests_inferred: 27324` — call-graph-reachable test edges (the 10× larger pool that makes diff-aware selection comprehensive, not just for hand-tagged tests).
- `verified: 376 / 444` equations have full theory→code→test chain.
- `orphan_code: 4035` (utility / non-math code, not a defect).

This means once the diff-base bug is fixed, `retest --scope=branch` should return a substantial `must_retest` set for our 179-commit diff. The selectivity (how much of the 2170 total can be safely skipped) is what we need to know to size the dev-workflow win.

## Related tools (also under-leveraged)

While filing the bug, propose using these tools more proactively in ORPHEUS workflow:

- `verification_audit` / `verification_gaps` / `verification_coverage` — V&V state queries.
- `provenance_chain` — equation → code → test trace.
- `trace_error` — failing test → suspect equations.
- `staleness` — doc drift detection (currently 0 stale, good — keep it).

ORPHEUS CLAUDE.md says "use Nexus first" but the actual workflow leans on pytest / grep / find. The retest bug fix is the unlock for promoting Nexus to default reach.

## Suggested upstream fix shape

Without source access to `sphinxcontrib-nexus`, the speculation:

1. Audit the git-diff invocation in `retest` and `detect_changes` for the `branch` scope. Likely a missing `--base` resolution or a hard-coded `HEAD`.
2. If hypothesis 3 holds (no stored pre-diff symbol set), the architectural fix is either (a) symbol-ID snapshot at each rebuild, or (b) two-graph diff. The latter is closer to "rebuild the graph at main, rebuild at HEAD, diff the symbol sets".

## Files for the next-session agent

- This memo: `.claude/agent-memory/numerics-investigator/nexus_retest_branch_scope_bug.md`
- Sphinxcontrib-nexus repo: `deOliveira-R/sphinxcontrib-nexus` (per user memory `project_graphify_fork.md`).
- ORPHEUS branch state: `refactor/sn-operator-algebra`, tip `01999e0`, 179 commits ahead of `main`.
