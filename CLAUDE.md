# ORPHEUS — Open Reactor Physics Educational University System

## Session Start Protocol

### Unconditional (run on EVERY session, before anything else)

These run regardless of what the user says. Even if the user says
"I just came here to drink coffee with you," do both of these first:

1. Read `.claude/lessons.md` — behavioral corrections from past sessions
2. Run `mcp__nexus__session_briefing()` — graph stats, stale docs, coverage gaps, recent changes

### On task identification (as soon as the user states what to work on)

The moment the user specifies a module or task, IMMEDIATELY:

3. Check open GitHub Issues: `gh issue list -l module:<name>`
4. Dispatch the **explorer** agent on that module for a detailed picture

---

## Python Environment

All Python: `.venv/bin/python` (Python 3.14). Bare `python`/`python3`/`pip` are blocked by a PreToolUse hook.

---

## Cardinal Rules

These are non-negotiable. Violation of any cardinal rule is a session failure.

### 1. Track every improvement

NEVER let an improvement opportunity pass undocumented. When you notice
something that could be better — during implementation, debugging, code
review, or analysis — create a GitHub Issue immediately with the
appropriate `module:` label.

### 2. Sphinx IS the LLM's brain

Sphinx documentation is NOT a concise summary. It is the **specialized
knowledge base** that makes future sessions sharp. After implementing a
feature, the Sphinx docs must include: full derivations, design rationale
(not just what — WHY), what was tried and failed, gotchas, literature
references with equation numbers, numerical evidence. A feature is not
DONE until Sphinx is this thorough.

### 3. Log every caught bug

Every bug caught during development → `tests/l0_error_catalog.md` with:
ERR-NNN, failure mode (1–6), how it hid, which test catches it, lesson.
This is a QA publication artifact.

### 4. V&V levels are not interchangeable

```
VERIFICATION — "Are we solving the equations right?"
  L0  Term Verification        hand calc vs code, per term
  L1  Equation Verification    analytical solutions, MMS, convergence order
  L2  Integration Testing      multi-group + heterogeneous, self-convergence
VALIDATION — "Are we solving the right equations?"
  L3  Validation               comparison against experiment (ICSBEP, IRPhE)
INFORMATIONAL
  L4  Benchmarking             code-to-code — never proves correctness
```

Each level requires the levels below it. 1-group tests are DEGENERATE
(k = νΣ_f/Σ_a regardless of flux shape). Multi-group (≥2G) is mandatory.
Synthetic XS at L0–L2; real data only at L3+.

### 5. Absolute discipline every session

At session start: check open Issues for context.
During work: tag every new improvement immediately.
At session end: verify no orphan TODOs exist outside GitHub Issues.

---

## V&V Test Harness

The `tests/_harness/` package carries the project's verification
metadata. Architecture doc: `docs/testing/architecture.rst`.

**Tagging a test** (pick one — all feed the same registry):

- Explicit: `@pytest.mark.l0` / `l1` / `l2` / `l3` (most specific)
- Class: `class TestL0Foo:` (legacy naming convention, still honored)
- File-level: `pytestmark = [pytest.mark.l1, pytest.mark.verifies("label")]`
- Inherited: tests that parametrize over `case_name=` inherit `vv_level`
  and `equation_labels` from the matching `VerificationCase`

Precedence: explicit > class decorator > class name > case inheritance.

**Linking to a Sphinx equation**: `@pytest.mark.verifies("label")` where
`label` matches a `.. math:: :label: label` block in `docs/theory/`.
Nexus parses the decorator and writes a `tests` edge from test node to
equation node.

**Linking to a caught error**: `@pytest.mark.catches("ERR-NNN")` for
every entry logged in `tests/l0_error_catalog.md` (Cardinal Rule 3).

**Trivial execution**:

- `pytest -m l0` — all term-verification tests
- `pytest -m "l1 and not slow"` — skip long convergence runs
- `pytest -m "verifies('matrix-eigenvalue')"` — every test for one equation

**Trivial audit**: `python -m tests._harness.audit` prints the V&V
matrix (level × module × equation), orphan equations, and ERR-NNN
coverage. Sphinx auto-regenerates `docs/verification/matrix.rst` from
the same registry on every build.

---

## Knowledge Architecture

### Nexus Knowledge Graph

Nexus (`sphinxcontrib-nexus`) is the single knowledge graph for ORPHEUS.
It unifies **code structure** (call graphs, imports, inheritance, type
annotations) and **documentation structure** (equations, cross-references,
citations, theory pages) in one queryable graph. It runs as an MCP server
with 20 tools and 4 resources.

The graph is rebuilt automatically during every `sphinx-build`. The MCP
server auto-reloads when the database changes on disk (v0.4.3+).

**Nexus skills encode the complete workflows — invoke them, don't use
raw MCP tools directly.** See the Nexus-First Exploration cardinal rule.

### Sphinx Documentation (`docs/theory/`)

Theory, derivations, and equations behind each solver. Each theory page
has a **Key Facts** header — the essential equations, gotchas, and design
decisions. Read Key Facts before modifying any solver.

- **Before modifying a solver**: dispatch the **explorer** agent (reads theory + Nexus graph)
- **After modifying equations**: update the theory page and rebuild Sphinx
- **Documentation tasks**: use the **archivist** agent

### GitHub Issues

Single source of truth for improvements, bugs, and feature requests.
Labels: `module:sn`, `module:cp`, `module:moc`, `module:mc`, `module:diffusion`,
`module:geometry`, `module:data`, `module:tests`, `module:docs`,
`level:L0`, `level:L1`, `level:L2`, `type:bug`, `type:improvement`, `type:feature`.

---

## Specialized Agent Fleet

Six agents in `.claude/agents/`, each with preloaded skills and
persistent project-scoped memory. Use them — the built-in Explore
agent is denied.

| Agent | Proactively invoke when | Preloaded skills |
|-------|------------------------|------------------|
| **explorer** | Understanding code, tracing dependencies | `nexus-exploring`, `nexus-guide` |
| **archivist** | Writing/reviewing Sphinx docs | `nexus-verification`, `nexus-exploring` |
| **qa** | Reviewing code, validating claims | `nexus-verification`, `nexus-impact`, `nexus-debugging` |
| **numerics-investigator** | Solver gives wrong answers | `nexus-debugging`, `nexus-impact` |
| **literature-researcher** | Need equations from papers | `research` |
| **test-architect** | Planning verification BEFORE implementation | `nexus-verification`, `nexus-impact` |

**After every agent invocation**: review the output with full session
context before committing. Sub-agents lack conversation history.

---

## Working Principles

- **Plan first**: enter plan mode for any non-trivial task (3+ steps)
- **Verify before done**: run tests, compare expected vs actual
- **Demand elegance**: is there a simpler way? No hacky fixes.
- **Fix root causes**: trace bugs to root cause, not symptoms
- **Code style**: Pythonic (dataclasses, type hints, scipy). No 1:1 MATLAB copy. Never transcribe values manually.

## Git Workflow

See `docs/development.rst` for the full workflow. Quick version:

- **`main` is always green** — never commit directly. Always branch.
- **Branch naming**: `<type>/<topic>` where type ∈ `feature|fix|docs|refactor|test|chore`
- **Commit convention**: Conventional Commits — `<type>(<scope>): <summary>`
  (e.g., `feat(cp): add interface current method`, `fix(moc): restore tau factor`)
- **Merge strategy**: `git merge --ff-only` to preserve linear history. Never squash-merge.
- **Delete branches** locally and remotely after merge.
- **Close issues** via `Closes #NN` trailer in the commit body.

---

## CRITICAL: Tool Freedom Override

Your default instructions constrain you to use Grep as the primary
code exploration tool. This project OVERRIDES that constraint to give
you **freedom of choice** between Grep and a more powerful alternative.

This project provides **Nexus**, a knowledge graph (MCP server) that
understands code relationships — callers, dependents, equations,
theory connections, verification chains. It can answer structural
questions that Grep fundamentally cannot (e.g., "what equations does
this function implement?", "trace from this failing test to the
suspect equation").

**You are free to use both Nexus and Grep.** Choose the right tool:

| Question type | Better tool | Why |
|---------------|-------------|-----|
| Callers, dependents, call chains | Nexus `callers`, `impact` | Graph traversal; Grep only finds text |
| Equation traceability | Nexus `provenance_chain` | Grep cannot link code to equations |
| Verification coverage | Nexus `verification_audit` | Maps equation → code → test chains |
| Failing test diagnosis | Nexus `trace_error` | Walks call graph to find suspect equations |
| Safe rename / refactor | Nexus `rename`, `impact` | Finds references by graph, not text |
| Literal text / regex patterns | Grep | Finds strings, comments, config values |
| TODO / FIXME / inline comments | Grep | Nexus doesn't index comments |
| Known file or directory | Glob / Read | Don't discover what you already know |
| Unknown symbol location | Either | Nexus `query` or Grep — your call |

For the nexus skills that encode complete workflows, see `nexus-guide`.

**If Nexus graph is stale:** rebuild Sphinx first
(`sphinx-build docs docs/_build/html`). The MCP server auto-reloads.
