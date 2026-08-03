# Delegation & sub-agent dispatch

Companion to CLAUDE.md **Cardinal Rule 5** (proactively delegate to sub-agents). Rule 5
establishes *that* and *why* you delegate; this rule covers the posture and the exceptions.
The dispatch *mechanics* (DISPATCH_REQUEST/RESULT bridging) live in the
`subagent-handoff-protocol` skill.

## Dispatch freely — don't ask permission

Dispatch sub-agents (`explorer`, `numerics-investigator`, `literature-researcher`, `qa`,
`archivist`, `test-architect`, …) freely during investigation and implementation. Do NOT
pause to ask "shall I dispatch X?" — just do it. Use parallel dispatches for independent
investigations. (The user values momentum and parallelism; the fleet's existence is already
approved.) Still **review every agent's output with full session context before committing.**

## Exception — surgical, high-correctness carves: main agent writes directly

For **surgical / high-importance** work — operator-algebra carves, convention changes
crossing ≥3 subsystems, anything in the `refactor/sn-operator-algebra` family — do **NOT**
dispatch `method-implementer`. The main agent writes the code directly with the user
steering step-by-step.

- Why: the user has since-inception codebase knowledge that corrects the implementation in
  real time; `method-implementer` runs a brief to completion, and for surgical work the loss
  of turn-by-turn correction costs more than the parallelism gains (user, 2026-05-20).
- The constraint is specifically on `method-implementer` (batch code production). Other
  agents (`test-architect`, `explorer`, `qa`, `archivist`) remain available and encouraged.
- Default for any "surgical" framing the user uses: direct main-agent writes + frequent
  `AskUserQuestion` checkpoints. When routine refactor cycles resume, `method-implementer`
  re-enters scope at the user's signal.

## Briefing a literature pull — check the local folder FIRST

When dispatching `literature-researcher` to acquire a paper, the brief MUST direct it to
check `scratch/literature/` (spell out the exact path) **before** any online search. The
user maintains this folder actively and has ALL Nuclear Science & Engineering volumes
locally.

- Prefer "extract equations from the local PDF at `<full path>`" over "find and acquire
  paper X" — the first is unambiguous and cheap.
- If the paper is NOT in the folder, the agent's first response must be "not in local
  folder; acquire it, or will you add it?" — NOT a unilateral pivot to a secondary source.
- Pivoting to a secondary source is a structural decision (different math path, possibly
  weaker verification claim) that needs **user approval**, not agent autonomy.

### OCR sidecars — the search/quote surface; the scan stays the SSOT

Every library PDF has a Mistral-OCR sidecar at `scratch/literature_ocr/<stem>.md`
(per-page `## p. <N>` sections, N = 1-based PDF page = the number `Read pages=` takes;
the provenance header carries the printed-page mapping). Regenerate/extend with
`.venv/bin/python tools/ocr_literature.py` (cache-idempotent; re-runs never re-bill).

- **Sidecar FIRST**: grep/read the sidecar before any visual page read — it is the
  search, navigation, and prose-quotation surface. Visual reads (`Read` with `pages=`)
  are for verification and for figures/layout the OCR can't carry.
- **The scan is the SSOT**: any **load-bearing equation** is spot-verified against the
  rendered page before it is transcribed into a theory page or a solver (the
  ERR-032-class hazard: a plausible-but-wrong reference equation).
- **Output discipline (the 2026-07-22 content-filter lesson)**: two agents died on
  `400 Output blocked by content filtering policy` while *generating* long verbatim
  transcriptions of scanned nuclear literature. Paraphrase + page-cite; keep verbatim
  quotes short; SELECT from the sidecar text instead of transcribing from page renders;
  and write findings to the deliverable file INCREMENTALLY so a mid-run kill preserves
  progress (both dead agents lost everything — zero bytes written).
- **Zotero liveness — 0 hits is not "empty".** A Zotero MCP server returning **0 hits on a
  known-present item** *together with* **connection-refused on port 23119** is BROKEN, not
  a library that lacks the paper. The misreading is silent and burns turns on repeated
  0-hit queries. Stop querying, record "Zotero down — no annotations checked", and proceed
  on the local folder + Tier-2. Applies to any agent briefed to consult Zotero, not just
  `literature-researcher`.
