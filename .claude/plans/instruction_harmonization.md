# Instruction-architecture harmonization — Session-1 record

**Status: DONE (Session 1, 2026-06-21).** Branch `chore/instruction-harmonization`,
ff-merged to main. This is the durable record + the hand-off to Session 2.

## What this was

The `chore/instruction-architecture` branch (4 commits, ~5 days stale vs 116 main
commits) introduced a new instruction *architecture* — `.claude/rules/*` (always-on
**floor**) + on-demand skills (**ceiling**) + a slimmed CLAUDE.md (247→177), tuned to
the Opus-4.8 system prompt. Rather than a mechanical `git rebase` (which would replay
stale diffs), it was **re-authored** onto current main: branch-as-spec, main-as-content-truth.

## THE harmonization heuristic (reusable)

Favor by **dimension**, not by file:
- **Directive / harness-architecture** → the directive-tuning tree wins (here: the branch). CLAUDE.md structure, `rules/`, skill *directives*, tool-routing.
- **Accumulated knowledge** (lessons, ERR catalog, vv failure-modes) → the knowledge-accumulating tree wins (here: main).
- **Neither clearly ahead** → judge per-item by relevance.
- **Cross-reference modernization is universal** — apply the graduated-link rewiring (`[[feedback-*]]` → `rules/`/skills) to every file regardless of which side won the content. Content-favoring and link-modernization are independent.

(Also saved as a main-agent memory; Session 2 may graduate it into a rule/skill.)

## What landed (5 commits on main)

1. `f144710` — instruction layer (12 files): slim CLAUDE.md + 5 new `rules/*` + nexus-elegance + coding-elegance/behavioral-auto-regression/nexus-exploring updates + lessons.md (main's L1–L28 + branch's link modernization) + tool_routing_ablation_study plan. CLAUDE.md relocation (Tool-Freedom→nexus-tools, V&V-Harness→vv-testing) **qa-verified LOSSLESS**.
2. `e9c9c3f` — agent-memory reconcile (90 files): committed the campaign closeouts/specs across all 9 agent dirs (triaged COMMIT-all by 3 archivist passes); fixed index coherence (+5 missing pointers, −9 dead orphans; 0 orphans remain).
3. `0de6a05` — V&V knowledge accumulators (vv-principles/SKILL.md + error_catalog.md).
4. `119422d` — SessionStart environment-health gate (session-health.sh + settings.json wiring + session-start.txt protocol).
5. (this commit) — durable plans: this record + the deferred pyright/operator-generics plans.

## Deferred — Session 2 (post-compact MEMORY.md deep-clean)

The user's Session-2 task: deep-clean the main-agent `MEMORY.md` (27.2 KB > 24.4 KB
limit) — harvest durable content into skills/rules, clean prose, bias toward recent
sessions, cut context noise. While there, also resolve these **flagged findings**:
- ~~**Pre-existing dangling "Cardinal Rule 6"** ref in vv-principles/SKILL.md + error_catalog.md (CLAUDE.md has only Rules 1–5). Decide the correct target (the ≥2G / 1-group-degeneracy directive) or drop the citation.~~ **RESOLVED (Session 2, user chose "retire the shorthand"):** "Cardinal Rule 6" never existed in CLAUDE.md (`git log -S` empty) — it was a long-dead shorthand for the 1-group-degeneracy rule, whose canonical home is now `vv-principles` SKILL.md §"1-group degeneracy — canonical statement" (anti-pattern #3 = operational form). Retired the shorthand across **21 files / 29 refs** (production tests, `docs/theory/*.rst`, `orpheus/derivations/`, agent-memory, archived plans) → "the 1-group-degeneracy rule". vv-principles anchor reworded to record the retirement + drop its false self-citation; one misattribution fixed (`discrete_ordinates.rst` "V&V harness wiring" was never the degeneracy rule). Sphinx clean; grep=0.
- **6 per-agent `AGENT.md` `## CRITICAL: Tool Freedom Override` sections** now overlap `rules/nexus-tools.md`. Retiring them needs the fact: *do `.claude/rules/` auto-load into sub-agent contexts?* If yes, the per-agent copies are redundant (twin-path, anti-pattern #1) → retire; if no, keep them as the sub-agent delivery mechanism.
- **Agent-index bloat**: method-implementer/MEMORY.md (+others) entries run far past the ~200-char one-line convention; qa/lessons.md has duplicate `L-NNN` heading numbers; numerics-investigator MEMORY.md lessons-pointer says `L1..L7` (now L8). The owning agents (or the clean pass) renumber/trim.

## Deferred — Session 3+ (pyright #226, unchanged)

The pyright signal-cleanup campaign: `.claude/plans/pyright_signal_cleanup.md` (authoritative).
Baseline at hand-off: **502** errors in `orpheus/` (260 production + 242 SymPy), 2353 full-repo.
The ratchet (`tests/test_pyright_ratchet.py`) is on main but baselined STALE at 706 → first
Session-3 step is `python -m tests._harness.pyright_ratchet --update`. Then the burn-down
(directory-by-directory, no scope exclusions, per-cluster commit + re-measure; B4
`LinearOperator[V]` generic is the highest-collapse keystone — design map at
`.claude/agent-memory/explorer/issue_226_operator_generics_map.md`). LSP-noise
investigation (option 2: explorer + a proposed `regen-pyrightconfig.sh` change for the
user to commit) also follows the merge.
