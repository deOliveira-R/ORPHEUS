#!/usr/bin/env bash
# PostToolUse hook for EnterWorktree — make a fresh worktree fully usable
# with zero manual protocol steps:
#
#   1. per-checkout venv + editable install (scripts/setup.sh, uv-first;
#      required for faithful pyright runs from the worktree root — see
#      memory `pyright-lsp-rooting` and ORPHEUS#226 context)
#   2. Sphinx/Nexus graph build inside the worktree (so the knowledge
#      graph is switchable the moment setup finishes)
#
# Registered with `asyncRewake` and ALWAYS exits 2: the agent is woken
# with this script's stderr as a system reminder — on success it should
# immediately call mcp__nexus__use_workspace(<worktree name>); on
# failure it sees the log tail and can react.
#
# Idempotent by delegation: scripts/setup.sh is safe to re-run, so
# re-entering an existing worktree just refreshes it (incremental
# sphinx build; venv untouched). Extra args are passed through to
# setup.sh (used by the standalone test: --no-docs).
set -uo pipefail

payload=$(cat)

# The hook JSON shape for EnterWorktree's tool_response is not pinned by
# docs — extract any absolute path under .claude/worktrees/ from the raw
# payload instead of guessing field names.
wt=$(WT_PAYLOAD="$payload" /usr/bin/python3 -c '
import os, re
m = re.search(r"(/[^\"\\\s]+/\.claude/worktrees/[^/\"\\\s]+)", os.environ["WT_PAYLOAD"])
print(m.group(1) if m else "")
')

if [ -z "$wt" ] || [ ! -d "$wt" ]; then
    # Fallback: the most recently created worktree of this project.
    base="${CLAUDE_PROJECT_DIR:-.}/.claude/worktrees"
    wt=$(ls -1td "$base"/*/ 2>/dev/null | head -1 | sed 's:/$::')
fi

if [ -z "$wt" ] || [ ! -d "$wt" ]; then
    echo "worktree-setup hook: could not locate the new worktree (payload had no .claude/worktrees path and $base is empty)." >&2
    exit 2
fi

# Keep pyright's executionEnvironments in sync with the worktree set
# (cross-tree type-identity noise — see regen-pyrightconfig.sh). Takes
# effect at the NEXT pyright-langserver start; cheap, so run it before
# the long venv/Sphinx build.
bash "$(dirname "$0")/regen-pyrightconfig.sh" || true

log="${TMPDIR:-/tmp}/orpheus-worktree-setup-$(basename "$wt").log"

if (cd "$wt" && bash scripts/setup.sh "$@") >"$log" 2>&1; then
    echo "Worktree ready: $wt — venv + editable install + Nexus graph built (log: $log). Switch the knowledge graph NOW with mcp__nexus__use_workspace(\"$(basename "$wt")\")." >&2
else
    echo "Worktree setup FAILED for $wt — tail of $log:" >&2
    tail -15 "$log" >&2
fi
exit 2
