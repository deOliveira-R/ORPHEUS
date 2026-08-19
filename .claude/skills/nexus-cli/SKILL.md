---
name: nexus-cli
description: "Use when the user needs to run Nexus from the shell — index a repo, start the MCP server, ask where the graph lives, or run a query/audit outside an agent session. Examples: \"Index this repo\", \"Analyze the codebase\", \"Start the nexus server\", \"Where is the graph?\", \"Run the V&V audit from CI\""
---

# Nexus CLI

Every MCP tool has a CLI twin, so anything an agent can ask, a shell script
or a CI job can ask too. `nexus --help` lists all 47 subcommands.

## Where the graph lives — DERIVE it, never declare it

⛔ **The database path is a CONVENTION, not a setting: `<project root>/.nexus/graph.db`.**
There is no `[graph].db` key — it was retired precisely so two surfaces
cannot disagree about where the graph is. The CLI, the MCP server and the
Sphinx extension each derive it from the project root.

```bash
nexus config db                     # the resolved path, bare, for shell capture
nexus config db --project-root DIR  # anchored to a checkout (a worktree, say)
nexus config                        # every resolved setting, as `key = value`
```

Exit 0 = the value is known; exit 1 = the key is unset or unknown — so
`nexus config db` is safe to use in a conditional.

**Keys `nexus config` will print:** `root`, `config`, `db`, `output`,
`scope.prefixes`, `catalog.errors`.

⟹ **Pass `--db` only to override the convention deliberately** — a scratch
graph, a second checkout, a graph you are diffing against. Retyping the
path by hand creates a second declaration, and the second declaration is
the one that goes stale. (It did: a hook with a hardcoded path kept a
`[ -f ]` guard that turned the drift into a quiet `exit 0`,
indistinguishable from "this project has no graph".)

The config file is found by walking UP from `--project-root` on the
subcommands that accept it, and from the **current working directory** on
those that do not. So run `nexus` from inside the tree, pass
`--project-root` where it is accepted, or pass `--db`.

## What `.nexus/config.toml` actually configures

Not the database. Five tables, and an unknown key is REPORTED rather than
silently dropped:

| table | keys |
|---|---|
| `[graph]` | `output`, `extra_source_dirs`, `exclude_patterns`, `analyze_tests`, `test_patterns`, `infer_implements`, `verification_registry`, `max_viz_nodes`, `git_timeout_seconds` |
| `[scope]` | `prefixes` |
| `[catalog]` | `errors` |
| `[replies]` | `max_characters`, `items_per_list`, `items_per_brief_line`, `neighbors_per_edge_type`, `nodes_per_impact_depth` |
| `[briefing]` | `project_hubs`, `stale_pages`, `symbols_per_stale_page`, `coverage_gaps` |

`[graph].output` is where the interactive `graph.html` is written, relative
to the Sphinx HTML output — that one artefact has to live in the HTML tree
because the docs link and iframe it. `graph.db` and `graph.json` do not.

⭐ **`[replies]` decides what a session can afford.** A reply lands in an
agent's context and stays there; these are the numbers to raise as context
windows grow.

## Runtime traces live beside the graph, and SURVIVE the build

`.nexus/traces/` — deliberately not under the Sphinx output. A profiled or
context-carrying suite run costs minutes to reproduce, so it is durable
state, while `graph.db` is rewritten on every build. Under `docs/_build/`
a single `rm -rf` destroyed it. **A directory's lifetime is its
most-derived member's.**

## Commands

```bash
nexus setup                    # install skills + rules; --check / --diff / --force
nexus analyze .                # AST-index a tree INTO the configured graph (merges)
nexus serve                    # MCP server on the configured graph
nexus status                   # graph summary
nexus briefing                 # full session briefing (JSON)
nexus workspaces               # every checkout of this project and its graph
```

`analyze` takes `--sys-path A B`, `--auto-sys-path` (numbered dirs),
`--exclude`, `--json PATH`, and `--db` for a scratch graph.

**Structure:** `query`, `context`, `neighbors`, `callers`, `callees`,
`impact`, `shortest-path`, `graph-query`, `processes`, `file-brief`,
`rename`.
**Smells:** `twin-paths`, `discriminations`, `native-place`,
`dead-functions`, `protocol-conformers`, `communities`, `bridges`,
`god-nodes`.
**Docs ↔ math:** `provenance`, `coverage`, `audit`, `gaps`, `doc-impact`,
`staleness`, `dead-references`, `errors`, `trace`.
**Runtime:** `runtime-ingest`, `runtime-runs`, `runtime-hotspots`,
`runtime-edges`, `runtime-branches`, `runtime-exercisers`,
`runtime-timeline`.
**Change:** `changes`, `retest`, `migration`, `ingest`, `visualize`.

### ⭐ `--run` — the flag that lets the CLI CONTRADICT a coverage claim

`retest`, `coverage` and `audit` each take `--run` (comma-separated runs
are unioned). Given one, every claim gains an execution verdict:
`corroborated` / `refuted` / `out_of_capture` / `no_implementation`.
Without one, nothing can contradict a claim at all.

```bash
nexus coverage --run geom_ctx,num_ctx --status verified
nexus audit --run geom_ctx,num_ctx --group-by module
nexus retest --scope branch --run geom_ctx
```

⛔ All three **refuse** a run that cannot carry attribution — a `cprofile`
run has no `exercised_by`, and answering from it would report "nothing is
covered" rather than "wrong instrument". The refusal names the runs that
DO qualify, and exits 1.

⚠ `exercised_by` needs the capture to have been taken with
`dynamic_context = test_function` **and** an unbroken `__init__.py` chain
from the rootdir to the test file — coverage names its context from the
MODULE, the graph names its node from the PATH, and a broken chain makes
them disagree silently. Check `nexus runtime-runs` for `exercised_by` in a
run's families before trusting any verdict.

### Exit codes for CI

`dead-references` and `errors` take `--exit-code` (non-zero when anything
is found) and `--quiet-when-clean`, so they drop straight into a pipeline.

## Sphinx integration

```python
extensions = ['sphinxcontrib.nexus']   # docs/conf.py
```

Settings belong in `.nexus/config.toml`, not `conf.py` — one declaration,
read by the extension, the CLI and the server alike. `conf.py` still
accepts the older `nexus_*` values as a LOWER-precedence tier, but settings
split across two files are settings that drift.

The graph rebuilds on every `sphinx-build`. For standalone use, re-run
`nexus analyze .` — and restart the MCP server, which loads from the
database at startup.
