#!/bin/bash
# Push dead documentation references into an agent's context.
#
# WHY A HOOK AND NOT A SKILL: steering an agent to go looking for dead
# references is probabilistic — measured, it depends on how the user
# happens to phrase the request, and an agent that is never asked will
# never look. A dead reference draws NO Sphinx warning at any severity,
# so if the agent doesn't look, nothing reports it. Pushing the finding
# is deterministic. Use this for the drift you cannot afford to miss.
#
# Wire it up in .claude/settings.json — SessionStart gives every session
# the current state; PostToolUse on Edit|Write catches drift you just
# created:
#
#   "hooks": {
#     "SessionStart": [
#       { "hooks": [ { "type": "command",
#                      "command": ".claude/hooks/nexus-dead-refs.sh" } ] }
#     ]
#   }
#
# Failure contract, same as nexus-brief.sh: quiet exit 0 on ANYTHING
# unexpected (no graph, no binary, corrupt db). Ambient information must
# never block or noise a session.

set -u

root="${CLAUDE_PROJECT_DIR:-$(pwd)}"
limit="${NEXUS_DEAD_REFS_LIMIT:-15}"

# The binary must be found BEFORE the graph, because it is what knows
# where the graph is. Asking `nexus config db` instead of hardcoding a
# path is not a style choice: the location is a convention derived from
# the project root (`<root>/.nexus/graph.db`), and a copy here is a
# second declaration that drifts silently. It did — when the store moved
# out of the Sphinx build output, this hook's hardcoded
# `docs/_build/html/_nexus/graph.db` stopped existing and the `[ -f ]`
# guard below turned the mistake into a quiet `exit 0`, indistinguishable
# from "this project has no graph".
nexus_bin="$root/.venv/bin/nexus"
[ -x "$nexus_bin" ] || nexus_bin="$(command -v nexus 2>/dev/null)"
[ -n "$nexus_bin" ] && [ -x "$nexus_bin" ] || exit 0

db="${NEXUS_DB:-$("$nexus_bin" config db --project-root "$root" 2>/dev/null)}"
[ -n "$db" ] || exit 0
[ -f "$db" ] || exit 0

# --quiet-when-clean: a clean project must cost zero context, or the
# hook trains agents to skim past it.
out=$("$nexus_bin" dead-references \
        --db "$db" \
        --format text \
        --limit "$limit" \
        --quiet-when-clean 2>/dev/null) || exit 0
[ -n "$out" ] || exit 0

printf '%s' "$out" | python3 -c '
import json, sys
print(json.dumps({
    "hookSpecificOutput": {
        "hookEventName": "SessionStart",
        "additionalContext": sys.stdin.read(),
    },
}))
' 2>/dev/null
exit 0
