#!/usr/bin/env python3
"""PreToolUse hook on Edit, Write and MultiEdit, set in the front matter of an agent whose role
writes nothing tracked (a reviewer, a Support agent): refuses a write outside the scratch tree, the
temporary directory and the agent's own memory. Usage, as the front matter's hook command:
``python3 .claude/hooks/write-scope.py <agent-name>``.

The brief template's "read-only" means no edit to a tracked file (the workflows page, "The brief");
for these agents the harness enforces it rather than the prose asking for it (instrument-doctrine X3).
The agent's own memory is always writable (the agent-definitions plan, ruling R2). Not seen, stated so
a pass is read for what it is: a write through Bash, which the agent's definition forbids and no hook
parses. Exit 2 blocks the call and hands the reason back; anything else lets it through.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path


def allowed_roots(project: Path, agent: str) -> list[Path]:
    temps = {Path(tempfile.gettempdir()), Path("/tmp"), Path("/private/tmp")}
    return [project / "scratch", project / ".claude" / "agent-memory" / agent, *(t.resolve() for t in temps)]


def main() -> int:
    agent = sys.argv[1] if len(sys.argv) > 1 else ""
    try:
        payload = json.load(sys.stdin)
    except Exception:
        return 0
    if payload.get("tool_name") not in {"Edit", "Write", "MultiEdit"} or not agent:
        return 0
    raw = (payload.get("tool_input") or {}).get("file_path") or ""
    if not raw:
        return 0
    project = Path(os.environ.get("CLAUDE_PROJECT_DIR") or payload.get("cwd") or ".").resolve()
    path = (project / raw).resolve() if not os.path.isabs(raw) else Path(raw).resolve()
    if any(path.is_relative_to(root) for root in allowed_roots(project, agent)):
        return 0
    print(f"write-scope: {agent} writes only under scratch/, the temporary directory and its own memory "
          f"(.claude/agent-memory/{agent}/); {path} is outside them. Put the content in your report or in a "
          f"scratch file, or ask the orchestrator with SendMessage.", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
