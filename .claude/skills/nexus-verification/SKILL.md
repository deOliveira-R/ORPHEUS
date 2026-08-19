---
name: nexus-verification
description: "Use when the user wants to check verification status, documentation coverage, or doc-code drift. Examples: \"What's verified?\", \"Which docs are stale?\", \"What equations have no tests?\", \"Do the docs reference things that no longer exist?\", \"Documentation coverage report\""
---

# Verification & Documentation Quality with Nexus

IMPORTANT: This skill is the dedicated tool for V&V assessment and
documentation drift detection.
Use Nexus for structural queries; use Grep freely for text search.

## Workflow

```
1. verification_audit()                                   → Coverage + drift, one call
2. verification_coverage({status_filter: "implemented"})  → Gaps: code but no tests
3. verification_gaps()                                    → Untagged tests, unverified equations
4. dead_references()                                       → Docs citing symbols that NO LONGER EXIST
5. staleness()                                             → Docs older than the code they document
```

## Two kinds of doc drift — they need different tools

- **Stale** (`staleness`): the doc still describes something real, but the
  code moved on. Git timestamps decide.
- **Dead** (`dead_references`): the doc references a class, function,
  attribute, or equation label that is **gone**. Sphinx renders a dead
  reference as plain text and emits NO warning at any severity, so nothing
  else in the toolchain will tell you. Report each dead target with the
  sites that still reference it; project-rooted names only, with re-export
  and inheritance rescue so live-but-indirect references aren't flagged.

## Coverage Status Values

| Status | Meaning | Action |
|--------|---------|--------|
| **verified** | Equation + code + test | Fully traced |
| **tested** | Code + test, no equation link | Add IMPLEMENTS doc |
| **implemented** | Equation + code, no test | Write a test |
| **documented** | Equation only, no code | Implement or mark future |
| **no_implementation** | An author DECLARED that nothing implements it | Nothing — this is the answer, not a gap |
| **orphan_code** | Code with no equation | Document the theory |

## ⭐⭐ A coverage claim can be REFUTED — pass `run=`

Every `verifies` marker is a **claim**, authored, stamped `confidence=1.0`,
pointing at an equation rather than at code. Nothing in the static graph can
contradict one. A coverage run joined to the graph can:

```
verification_audit({run: "geom_ctx,num_ctx"})       # comma-separated = union
verification_coverage({run: "…", status_filter: "verified"})
retest({scope: "branch", run: "…"})
```

Each claim then carries an `execution` verdict:

| verdict | meaning |
|---|---|
| `corroborated` | the claiming test EXECUTED an implementing node |
| `refuted` | the claimant ran in this capture and touched none |
| `out_of_capture` | the claimant is in no capture — **not evidence either way** |
| `no_implementation` | the equation has no implementing code to adjudicate against |

⛔ **The last two are not findings.** They say the claim could not be checked,
for two causes needing OPPOSITE repairs — a wider capture, and a declared
`implements` link. Report them separately or an audit reads as far worse (or
far better) than it is.

⛔ **A refutation is only as good as the CODE side.** Read `code_evidence` on
the row: `declared` means a directive or registry linked that code to that
equation; `inferred` means nobody did and the link is a shared-name-token
guess. `[M]` 2026-08-19 on one corpus, **13 121 of 13 508** `implements` edges
are still guesses — so most refutations are refuting a guess.

⛔ **`summary.claims_*` sums authored claims with two heuristic BFS tiers.**
`[M]` it printed 11 034 corroborated where the authored figure was 282. Split
on `TestReference.source == "declared"`, and note the aggregate INFLATES as
declarations land.

⚠ **A capture that binds nothing looks exactly like a capture that found
nothing.** `exercised_by` needs `dynamic_context = test_function` in the
coverage config AND an unbroken `__init__.py` chain from the rootdir to the
test file — coverage names the context from the MODULE, the graph names the
node from the PATH, and a broken chain makes them disagree silently. `[M]`
2026-08-19: one missing `__init__.py` made a 3330-test tree unattributable,
with `rc=0` and a plausible bind count. **Check `runtime_runs` for
`exercised_by` in the run's `families` before trusting any verdict.**

## Test Inventory Queries

Tests are indexed in the graph. Test node IDs: `py:function:tests.<file>.<function>`.

```
# List all test modules
query({text: "tests.test_", node_types: "module", limit: 50})

# List test functions in a module
neighbors({node_id: "py:module:tests.test_sn_1d", direction: "out", edge_types: "contains"})

# Find tests covering a function
impact({target: "py:function:orpheus.sn.solver.solve_sn", direction: "upstream"})

# Trace test → equations
trace_error({test_node_id: "py:function:tests.test_cp_slab.test_slab_cp_eigenvalue"})
```

## Checklist

- [ ] `verification_coverage` for full V&V status
- [ ] Review "implemented" entries (verification gaps)
- [ ] Review "documented" entries (unimplemented equations)
- [ ] Review "orphan_code" entries (undocumented theory)
- [ ] `dead_references` for docs citing symbols/equations that no longer exist
- [ ] `staleness` to find docs needing updates
- [ ] If a capture exists, re-run the audit with `run=` — a claim nothing
      executed is a gap the static view cannot see
- [ ] Create GitHub Issues for each gap found

See [../nexus-exploring/reference.md](../nexus-exploring/reference.md) for full reference.
