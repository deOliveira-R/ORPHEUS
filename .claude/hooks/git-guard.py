#!/usr/bin/env python3
"""PreToolUse hook on Bash: refuses the two git commands this repository forbids outright
(CLAUDE.md § How a session runs; the process-discipline rule).

1. ``git add -A`` / ``--all`` / a bare ``.`` without ``-u``: scratch/ is untracked and large;
   stage by path.
2. ``git commit`` while HEAD is ``main``: main is always green and receives only ff-merges.

Exit 2 blocks the call and hands the reason back to the agent; anything else lets it
through. ``git`` must sit at command position (line start, or after ; && || | $( ), so a
quoted mention such as ``grep "git add -A"`` is text, not a command, and passes. A compound
that branches and commits in one line is refused as well: branch in its own command first.
The payload arrives on stdin as JSON (``tool_name``, ``tool_input.command``, ``cwd``).
"""
import json
import re
import subprocess
import sys

AT_COMMAND = r"(?:^|[;&|(]|\n)\s*"


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except Exception:
        return 0
    if payload.get("tool_name") != "Bash":
        return 0
    cmd = (payload.get("tool_input") or {}).get("command") or ""
    cwd = payload.get("cwd") or "."
    for m in re.finditer(AT_COMMAND + r"git\s+add\s+([^;&|\n]*)", cmd):
        tokens = set(m.group(1).split())
        sweeps = {"-A", "--all", "--no-ignore-removal"} & tokens
        bare_dot = ({".", "./"} & tokens) and not ({"-u", "--update"} & tokens)
        if sweeps or bare_dot:
            print("git-guard: `git add -A|--all|.` is refused in this repository: scratch/ is untracked "
                  "and large. Stage by explicit path, or `git add -u` for tracked files, after reading "
                  "`git status --porcelain | grep -v '^?? scratch/'`.", file=sys.stderr)
            return 2
    m = re.search(AT_COMMAND + r"git\s+(?:-C\s+(\S+)\s+)?commit(?:\s|$)", cmd)
    if m:
        repo = m.group(1) or cwd
        try:
            head = subprocess.run(["git", "-C", repo, "symbolic-ref", "--short", "-q", "HEAD"],
                                  capture_output=True, text=True, timeout=5).stdout.strip()
        except Exception:
            head = ""
        if head == "main":
            print(f"git-guard: `git commit` while HEAD is `main` (in {repo}) is refused: main is always "
                  "green and receives only `git merge --ff-only`. Create a branch `<type>/<topic>` in "
                  "its own command first, then commit.", file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
