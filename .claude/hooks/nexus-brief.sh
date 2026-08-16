#!/bin/bash
# Edit-time Nexus brief — the graph's ambient push channel (D2).
#
# After an Edit/Write on a .py file, inject <=6 lines of graph context
# (blast radius, implemented equations, verifying tests, owning doc
# pages, staleness) into the agent's context via PostToolUse
# additionalContext — graph knowledge arrives WITH the edit, the way
# pyright pushes diagnostics, instead of waiting to be asked.
#
# Failure contract: quiet exit 0 on ANYTHING unexpected (file not in
# graph, no database, no venv, JSON oddities). The brief is ambient
# information; it must never block or noise an edit.

payload=$(cat)

parsed=$(printf '%s' "$payload" | python3 -c '
import json, sys
try:
    d = json.load(sys.stdin)
    print(d.get("tool_input", {}).get("file_path", ""))
    print(d.get("session_id", "anon"))
except Exception:
    pass
' 2>/dev/null) || exit 0
file=$(printf '%s\n' "$parsed" | sed -n 1p)
session=$(printf '%s\n' "$parsed" | sed -n 2p)

[ -n "$file" ] || exit 0
case "$file" in
  *.py) ;;
  *) exit 0 ;;
esac
case "$file" in
  /*) ;;
  *) file="$PWD/$file" ;;
esac

# The file's OWN checkout answers: walk up to the nearest checkout
# boundary (a .git entry — directory for a main tree, file for a
# linked worktree) and require the graph THERE. Never walk past it:
# a worktree file must not brief from the main checkout's graph
# (the L22 wrong-tree hazard) — no graph in this checkout means no
# brief, not a plausible-but-wrong one.
root=$(dirname "$file")
while [ "$root" != "/" ]; do
  [ -e "$root/.git" ] && break
  root=$(dirname "$root")
done
[ "$root" != "/" ] || exit 0
nexus_bin="$root/.venv/bin/nexus"
[ -x "$nexus_bin" ] || exit 0

# Ask nexus where this checkout's graph is rather than hardcoding it —
# the path is declared in .nexus/config.toml, and a copy here is a
# second declaration free to drift. `--project-root "$root"` keeps the
# worktree guarantee above intact: the answer is anchored to the
# checkout we just walked up to, never to the caller's directory.
db=$("$nexus_bin" config db --project-root "$root" 2>/dev/null) || exit 0
[ -f "$db" ] || exit 0

# Once per file per session: the first edit briefs, later ones stay
# quiet — repeated context injection is noise, not ambience.
dedup_dir="${TMPDIR:-/tmp}/nexus-brief-${session:-anon}"
mkdir -p "$dedup_dir" 2>/dev/null || exit 0
stamp="$dedup_dir/$(printf '%s' "$file" | /usr/bin/shasum | cut -d' ' -f1)"
[ -e "$stamp" ] && exit 0
: > "$stamp"

brief_out=$("$nexus_bin" file-brief "$file" \
  --db "$db" \
  --project-root "$root" 2>/dev/null) || exit 0
[ -n "$brief_out" ] || exit 0

printf '%s' "$brief_out" | python3 -c '
import json, sys
print(json.dumps({
    "hookSpecificOutput": {
        "hookEventName": "PostToolUse",
        "additionalContext": sys.stdin.read(),
    },
}))
' 2>/dev/null
exit 0
