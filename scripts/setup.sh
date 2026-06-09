#!/usr/bin/env bash
# One-shot dev-environment setup for ANY ORPHEUS checkout — the main clone, a
# git worktree, or a fresh cloud VM. Idempotent; safe to re-run.
#
#   ./scripts/setup.sh             # venv + editable dev install + Nexus graph
#   ./scripts/setup.sh --no-docs   # skip the (slow) Sphinx / Nexus-graph build
#
# WHY A VENV PER CHECKOUT (not a shared/symlinked one): `orpheus` is an EDITABLE
# install, so the venv must be created FROM this checkout — only then does
# `import orpheus` resolve to THIS branch's code. A venv symlinked from another
# checkout would silently import that other checkout's source instead.
#
# Cloud (Claude Code on the web): set the environment's "Setup script" field to
#   bash scripts/setup.sh
# It satisfies all three things the nexus MCP server needs in a fresh VM:
# a .venv, the installed `nexus` console script, and docs/_build/html/_nexus/graph.db.
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root of THIS checkout (works from a worktree too)

PY_VERSION=3.14           # matches the project; uv fetches it if absent

echo "==> Creating .venv + installing orpheus[dev] (editable, from $(pwd))"
if command -v uv >/dev/null 2>&1; then
    # uv hardlinks deps from a global cache → fast + ~zero extra disk per worktree.
    uv venv --python "$PY_VERSION"
    uv pip install -e ".[dev]"
else
    # Stdlib fallback (slower, no cross-venv dedup). Needs python >=3.11 on PATH.
    python3 -m venv .venv
    .venv/bin/python -m pip install --upgrade pip wheel
    .venv/bin/pip install -e ".[dev]"
fi

if [[ "${1:-}" != "--no-docs" ]]; then
    echo "==> Building Sphinx docs → writes the Nexus graph (docs/_build/html/_nexus/graph.db)"
    .venv/bin/python -m orpheus.derivations.generate_rst
    .venv/bin/python -m sphinx -b html docs docs/_build/html
fi

echo "✓ Done. import orpheus → this checkout. Activate with:  source .venv/bin/activate"
