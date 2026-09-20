#!/usr/bin/env bash
# InstructionsLoaded hook — appends one line per instruction file the harness injects
# (CLAUDE.md, .claude/rules/*.md; reason: session_start, path_glob_match, include, compact, ...).
# The payload is logged whole, so the log carries whatever fields the harness sends.
# Acceptance instrument for the harness restructure (.claude/plans/harness_context_budget.md, T2/T4).
root="${CLAUDE_PROJECT_DIR:-.}"
log="$root/scratch/_harness_eval/instructions_loaded.log"
mkdir -p "$(dirname "$log")"
{ printf '%s ' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"; tr -d '\n'; printf '\n'; } >> "$log"
exit 0
