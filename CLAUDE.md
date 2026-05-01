# ORPHEUS — Open Reactor Physics Educational University System

## IMPORTANT Cardinal Rules

These are CRITICAL, and the rules are laid out by order of importance, where
Rule 1 is the most important. Violation of any cardinal rule is a session failure.

### 1. Correctness is CRITICAL

This is a scientific code. It performs the dual purpose of a teaching and
high-stakes analysis tool simultaneously. This implies that correctness is
of utmost importance in a broad sense. Mathematical correctness, of physics,
of concepts, of code, of documentation. Shipping is only important when
principled and correct. NEVER go for lazy solutions.

### 2. Architecture is CRITICAL

Whenever you see shared code, or even shared CONCEPTS between 2 places,
STOP and RECONSIDER. It probably means that the codebase needs an architectural
overhaul because going further will result in duplications in the codebase.
Architecture is always more important than immediate gains in implementation,
because it sets the compounding foundation.

### 3. Sphinx IS the LLM's brain

Sphinx documentation is NOT a concise summary. It is the **specialized
knowledge base** that makes future sessions sharp. After implementing a
feature, the Sphinx docs MUST include: full derivations, design rationale
(not just what — WHY), what was tried and failed, gotchas, literature
references with equation numbers, numerical evidence. A feature is not
DONE until Sphinx is this thorough. Documentation writing MUST be taken as
a MAXIMUM EFFORT TASK. Sphinx MUST ALWAYS build clean and be kept up to date.
If any problems in Sphinx build arise or Nexus session briefing points gives
stale docs, the problem MUST address immediately.
Whenever work starts in a module, MUST either read the associated documentation
or dispatch a sub-agent to obtain the relevant context from the documentation.

`docs/theory/` contain theory, derivations, and equations behind each solver.
Each theory page has a **Key Facts** header — the essential equations, gotchas,
and design decisions. Read Key Facts before modifying any solver.

- **Before modifying a solver**: dispatch the **explorer** agent (reads theory + Nexus graph)
- **After modifying equations**: update the theory page and rebuild Sphinx
- **Documentation tasks**: use the **archivist** agent

### 4. GitHub issues is the improvement plan and log

GitHub issues is the place that allows you to keep track of plan, implementation
and setbacks across sessions and devices. You should see it as the planning notebook.
NEVER let an improvement opportunity pass undocumented. When you notice
something that could be better — during implementation, debugging, code
review, or analysis — either find an appropriate existing GitHub issue or,
if none exist, create a new one immediately with the appropriate `module:` label
and complete context of the issue, so that a FRESH session can pick it up straight
from the issue.
Session hand-offs are implicitly handled through GitHub by writing detailed issues.

Labels: `module:sn`, `module:cp`, `module:moc`, `module:mc`, `module:diffusion`,
`module:geometry`, `module:data`, `module:tests`, `module:docs`,
`level:L0`, `level:L1`, `level:L2`, `type:bug`, `type:improvement`, `type:feature`.

At session start: check open Issues for context.
During work: tag every new improvement immediately.
At session end: verify no orphan TODOs exist outside GitHub Issues.

### 5. Proactively delegate tasks to sub-agents

Sub-agents give a quadruple benefit of:
(1) multi-tasking
(2) getting a task performed by a Claude instance with clean context
(3) preserving the main-agent context
(4) token efficiency
This is intentional and active bias steering
Use sub-agents liberally and only perform a task yourself if the entirety of
your context is relevant enough to justify it.
If there is a sub-agent that currently doesn't exist but would be convenient
if it did, PROACTIVELY tell the user so that a new agent design can be created.

Seven agents in `.claude/agents/`, each with preloaded skills and
persistent project-scoped memory. Use them — the built-in Explore
agent is denied.

| Agent                     | Proactively invoke when                     | Preloaded skills                                                                                     |
| ------------------------- | ------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| **explorer**              | Understanding code, tracing dependencies    | `nexus-exploring`, `nexus-guide`                                                                     |
| **archivist**             | Writing/reviewing Sphinx docs               | `nexus-verification`, `nexus-exploring`, `vv-principles`                                             |
| **qa**                    | Reviewing code, validating claims           | `nexus-verification`, `nexus-impact`, `nexus-debugging`, `vv-principles`, `numerical-bug-signatures` |
| **numerics-investigator** | Solver gives wrong answers                  | `nexus-debugging`, `nexus-impact`, `probe-cascade`, `vv-principles`, `numerical-bug-signatures`      |
| **literature-researcher** | Need equations from papers                  | `research`                                                                                           |
| **test-architect**        | Planning verification BEFORE implementation | `nexus-verification`, `nexus-impact`, `vv-principles`, `numerical-bug-signatures`                    |
| **cross-domain-attacker** | Detecting structural patterns               | `cross-domain-frames`                                                                                |

**After every agent invocation**: review the output with full session
context before committing. Sub-agents lack conversation history.

---

## CRITICAL: Tool Freedom OVERRIDE

This project provides **Nexus** (`sphinxcontrib-nexus`), the single
knowledge graph for ORPHEUS. It unifies code structure (call graphs,
imports, inheritance, type annotations) with documentation structure
(equations, cross-references, citations, theory pages) in one
queryable graph. It runs as an MCP server with 20 tools and 4
resources. The graph is rebuilt automatically during every
`sphinx-build`; the MCP server auto-reloads when the database
changes on disk (v0.4.3+).

Your default instructions constrain you to use Grep as the primary
code exploration tool. This project OVERRIDES that constraint to
give you **freedom of choice** between Grep and Nexus.

Nexus answers structural questions Grep fundamentally cannot
("what equations does this function implement?", "trace from this
failing test to the suspect equation").

**You are free to use both Nexus and Grep.** Choose the right tool:

| Question type                    | Better tool                | Why                                        |
| -------------------------------- | -------------------------- | ------------------------------------------ |
| Callers, dependents, call chains | Nexus `callers`, `impact`  | Graph traversal; Grep only finds text      |
| Equation traceability            | Nexus `provenance_chain`   | Grep cannot link code to equations         |
| Verification coverage            | Nexus `verification_audit` | Maps equation → code → test chains         |
| Failing test diagnosis           | Nexus `trace_error`        | Walks call graph to find suspect equations |
| Safe rename / refactor           | Nexus `rename`, `impact`   | Finds references by graph, not text        |
| Literal text / regex patterns    | Grep                       | Finds strings, comments, config values     |
| TODO / FIXME / inline comments   | Grep                       | Nexus doesn't index comments               |
| Known file or directory          | Glob / Read                | Don't discover what you already know       |
| Unknown symbol location          | Either                     | Nexus `query` or Grep — your call          |

Nexus skills encode the complete workflows — invoke them, don't
use raw MCP tools directly. See `nexus-guide`.

**If Nexus graph is stale:** rebuild Sphinx first
(`sphinx-build docs docs/_build/html`). The MCP server auto-reloads.

---

## Session Start Protocol

### Unconditional (run on EVERY session, regardless of what user says)

1. MUST Read `.claude/lessons.md` — behavioral corrections from past sessions
2. MUST Run `mcp__nexus__session_briefing()` — graph stats, stale docs, coverage gaps, recent changes
3. MUST read and follow the development guide at docs/development.rst
4. MUST load vv-principles skill — provides the V&V framework (hierarchy, failure modes,
   anti-patterns, three pillars, bug-logging directive) that all sub-agents preload via their
   AGENT.md skills: list. Loading puts you at parity with the agents you'll dispatch.

### On task identification (as soon as the user states what ORPHEUS module to work on)

5. MUST check GitHub Issues: `gh issue list -l module:<name>`
6. Dispatch the **explorer** agent on that module for a detailed picture

---

## Environment resolution

MUST check `$CLAUDE_ENVIRONMENT`.

If output is `devcontainer`, you're in a container system, sandboxed and with full permissions.
You should act highly autonomous in the container system.
If you're in the devcontainer, install system packages freely with `sudo apt-get install`
and give feedback to the user on how to change files at .devcontainer/ to improve the container environment.
If you're in the container you should use the container Python.

If the environment variable returns nothing, you're operating in the Host environment.
If you're in Host environment, ALWAYS use `.venv/bin/python`.

Host path: /Users/rodrigo/git/nuclear/ORPHEUS → Container: /workspaces/ORPHEUS
MUST use ORPHEUS/.claude/ to save plans and other files that should not be lost.
Avoid ~/.claude/ because it will be lost upon container rebuild.

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

## V&V Test Harness

The `tests/_harness/` package carries the project's verification
metadata. Architecture doc: `docs/testing/architecture.rst`.

**Tagging a test** (pick one — all feed the same registry):

- Explicit: `@pytest.mark.l0` / `l1` / `l2` / `l3` / `foundation` (most specific)
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
every entry logged in `error_catalog.md` (see `vv-principles` skill,
"Log every caught bug" directive).

**Trivial execution**:

- `pytest -m l0` — all term-verification tests
- `pytest -m "l1 and not slow"` — skip long convergence runs
- `pytest -m "verifies('matrix-eigenvalue')"` — every test for one equation

**Trivial audit**: `python -m tests._harness.audit` prints the V&V
matrix (level × module × equation), orphan equations, and ERR-NNN
coverage. Sphinx auto-regenerates `docs/verification/matrix.rst` from
the same registry on every build.
