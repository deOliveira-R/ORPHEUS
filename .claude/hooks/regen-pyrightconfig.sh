#!/usr/bin/env bash
# regen-pyrightconfig.sh — keep pyright's view of git worktrees coherent.
#
# Problem (ORPHEUS#226 context / memory `pyright-lsp-rooting`): the
# pyright LSP plugin roots at the MAIN checkout for the whole session,
# even after EnterWorktree. A worktree file then resolves
# `from .geometry import X` inside the worktree but `from orpheus.…
# import Y` from the MAIN tree, so the same class acquires two type
# identities (".claude.worktrees.<name>.orpheus.….X" vs
# "orpheus.….X") and every call crossing them reports a false
# "X is not assignable to X". This is the known pyright nested-project
# dual-identity class (mypy errors on it as "Source file found twice
# under different module names"); the maintainer-sanctioned mechanism
# is execution environments (microsoft/pyright#10498).
#
# Fix: generate a machine-local (gitignored) pyrightconfig.json at the
# main root that mirrors pyproject.toml's [tool.pyright] and adds one
# executionEnvironment per existing worktree, rooting the worktree's
# subtree at ITSELF. Every import inside a worktree file then resolves
# to the worktree's own tree: the false dual identities disappear while
# every real diagnostic stays (verified A/B on sn/solver.py:
# 53 diagnostics / 18 false  →  35 diagnostics / 0 false).
#
# Caveat: pyright-langserver reads the config at STARTUP only, and the
# plugin client cannot deliver config-change notifications
# (anthropics/claude-code#27220, claude-plugins-official#1359) — a
# worktree created mid-session keeps the noise until the next session.
# Hence regeneration from three hooks (SessionStart, EnterWorktree via
# worktree-setup.sh, ExitWorktree): the config is always accurate for
# the NEXT server start, whichever event changed the worktree set.
#
# With no worktrees present the file is REMOVED so pyproject.toml's
# [tool.pyright] stays the single config source.
#
# Failure contract: quiet exit 0 on anything unexpected — this is
# ambient tooling and must never block a session or a hook chain.
set -uo pipefail

# Self-root via git so the script works from any cwd in any checkout
# (main, a worktree, or a hook's inherited directory).
common_dir=$(git rev-parse --path-format=absolute --git-common-dir 2>/dev/null) || exit 0
main_root=$(dirname "$common_dir")
[ -d "$main_root" ] || exit 0

py="$main_root/.venv/bin/python"
[ -x "$py" ] || py=python3

MAIN_ROOT="$main_root" "$py" - <<'EOF' 2>/dev/null
import json
import os
import pathlib
import sys

try:
    import tomllib
except ModuleNotFoundError:  # python < 3.11: skip silently
    sys.exit(0)

root = pathlib.Path(os.environ["MAIN_ROOT"])
cfg_path = root / "pyrightconfig.json"

worktrees = sorted(
    p for p in (root / ".claude" / "worktrees").glob("*")
    if (p / ".git").exists()
)
if not worktrees:
    cfg_path.unlink(missing_ok=True)
    sys.exit(0)

base = {}
pyproject = root / "pyproject.toml"
if pyproject.exists():
    base = tomllib.loads(pyproject.read_text()).get("tool", {}).get("pyright", {})

base["executionEnvironments"] = [
    {
        "root": str(p.relative_to(root)),
        "extraPaths": [str(p.relative_to(root))],
    }
    for p in worktrees
]
cfg_path.write_text(json.dumps(base, indent=4) + "\n")
EOF
exit 0
