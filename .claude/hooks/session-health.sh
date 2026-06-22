#!/usr/bin/env bash
# session-health.sh — SessionStart environment gate (#226 + Nexus readiness).
#
# Emits an "ENVIRONMENT HEALTH" block into the session context so the main
# agent can confirm a NOISE-FREE workspace BEFORE the session-start protocol
# runs — in particular BEFORE mcp__nexus__session_briefing(). Two gates:
#
#   PYRIGHT — the CLI must resolve first-party imports (proving the committed
#     pyrightconfig.json is correct, so the type checker is SIGNAL not noise).
#     The IDE/LSP diagnostic stream is ADVISORY; the CLI is the oracle (#226).
#     A canary smoke test on ONE leaf file: any reportMissingImports there means
#     the config is broken (first-party orpheus.* not resolving).
#   NEXUS — the graph DB must EXIST and be a readable sqlite with content, and
#     the nexus package must be installed (.venv/bin/nexus), BEFORE the agent
#     calls session_briefing (which would otherwise answer from a missing graph
#     or a dead server).
#
# For each gate: OK / PROBLEM + remediation / UNKNOWN (could-not-verify — never
# a false OK). The agent ACTS on it per session-start.txt: fix it (telling the
# user any manual step), or explain the problem and ask whether to proceed.
#
# Failure contract: NEVER block the session — bounded, quiet, always exit 0
# (the settings.json hook timeout bounds wall-clock). A check that cannot run
# prints UNKNOWN, not OK.

set -uo pipefail

root=$(git rev-parse --show-toplevel 2>/dev/null) || root="$PWD"
[ -d "$root" ] || root="$PWD"
py="$root/.venv/bin/python"
[ -x "$py" ] || py=python3

echo "=== ENVIRONMENT HEALTH (SessionStart gate — #226 pyright + Nexus) ==="

# ── Gate 1: PYRIGHT (the CLI is the oracle, not the LSP stream) ──────────────
canary="$root/orpheus/numerics/moment_layout.py"   # leaf, fast, first-party imports
pyright_status="UNKNOWN — npx/canary unavailable; pyright UNTRUSTED until a CLI run confirms clean."
if command -v npx >/dev/null 2>&1 && [ -f "$canary" ]; then
  out=$(cd "$root" && npx --no-install pyright --outputjson "$canary" 2>/dev/null) || out=""
  if [ -n "$out" ]; then
    missing=$(printf '%s' "$out" | "$py" -c 'import json,sys
try:
    d=json.load(sys.stdin)
    print(sum(1 for x in d.get("generalDiagnostics",[]) if x.get("rule")=="reportMissingImports"))
except Exception:
    print("err")' 2>/dev/null)
    case "$missing" in
      0)   pyright_status="OK — CLI resolves first-party imports (pyrightconfig.json correct). TRUST the CLI (npx pyright), NOT the streamed <new-diagnostics>: the IDE/LSP stream's import/'not defined'/cross-tree noise is a langserver artifact (#226), advisory only." ;;
      err|"") pyright_status="UNKNOWN — could not parse pyright output; verify by hand: npx pyright orpheus/numerics/moment_layout.py" ;;
      *)   pyright_status="PROBLEM — CLI reports ${missing} reportMissingImports on the canary → the config is NOT resolving first-party orpheus.*. REMEDIATION: (1) bash .claude/hooks/regen-pyrightconfig.sh ; (2) ensure the venv exists: bash scripts/setup.sh --no-docs ; (3) re-run the canary. Until clean, EVERY pyright import diagnostic is suspect." ;;
    esac
  else
    pyright_status="UNKNOWN — pyright did not run (npx cache/network?). Verify: npx --no-install pyright --version"
  fi
fi
echo "PYRIGHT: $pyright_status"

# ── Gate 2: NEXUS (loaded + executable BEFORE session_briefing) ─────────────
db="$root/docs/_build/html/_nexus/graph.db"
nexus_bin="$root/.venv/bin/nexus"
nexus_status="UNKNOWN"
if [ ! -f "$db" ]; then
  nexus_status="PROBLEM — no graph DB at docs/_build/html/_nexus/graph.db. session_briefing will answer from nothing. REMEDIATION: sphinx-build docs docs/_build/html (the build writes the graph; the MCP server auto-reloads)."
elif [ ! -x "$nexus_bin" ]; then
  nexus_status="PROBLEM — graph DB exists but the nexus package is not installed (.venv/bin/nexus missing). REMEDIATION: bash scripts/setup.sh --no-docs (installs sphinxcontrib-nexus into .venv)."
else
  db_ok=$("$py" -c 'import sqlite3,sys
try:
    c=sqlite3.connect("file:%s?mode=ro"%sys.argv[1], uri=True)
    n=c.execute("select count(*) from sqlite_master").fetchone()[0]
    print("ok" if n>0 else "empty")
except Exception:
    print("bad")' "$db" 2>/dev/null)
  case "$db_ok" in
    ok)    nexus_status="OK — graph DB present + readable, nexus package installed. session_briefing is safe to call. If mcp__nexus__* surface as deferred, ONE ToolSearch(select:mcp__nexus__session_briefing,...) loads them (deferral is NOT unavailability). Staleness/branch-mismatch, if any, is reported INSIDE the briefing's workspace block." ;;
    empty) nexus_status="PROBLEM — graph DB is empty (no tables). REMEDIATION: sphinx-build docs docs/_build/html to (re)build the graph." ;;
    *)     nexus_status="PROBLEM — graph DB exists but is unreadable/corrupt. REMEDIATION: rebuild via sphinx-build docs docs/_build/html; if that fails, bash scripts/setup.sh --no-docs." ;;
  esac
fi
echo "NEXUS: $nexus_status"
echo "=== end ENVIRONMENT HEALTH ==="
exit 0
