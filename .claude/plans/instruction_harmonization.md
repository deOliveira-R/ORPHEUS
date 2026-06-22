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

## Session 2 — DONE (2026-06-22)

Deep-cleaned the WHOLE memory substrate (not just main MEMORY.md). 17 commits,
ff-merged to main. The 3 flagged findings, all RESOLVED:
- ~~**Pre-existing dangling "Cardinal Rule 6"** ref in vv-principles/SKILL.md + error_catalog.md (CLAUDE.md has only Rules 1–5). Decide the correct target (the ≥2G / 1-group-degeneracy directive) or drop the citation.~~ **RESOLVED (Session 2, user chose "retire the shorthand"):** "Cardinal Rule 6" never existed in CLAUDE.md (`git log -S` empty) — it was a long-dead shorthand for the 1-group-degeneracy rule, whose canonical home is now `vv-principles` SKILL.md §"1-group degeneracy — canonical statement" (anti-pattern #3 = operational form). Retired the shorthand across **21 files / 29 refs** (production tests, `docs/theory/*.rst`, `orpheus/derivations/`, agent-memory, archived plans) → "the 1-group-degeneracy rule". vv-principles anchor reworded to record the retirement + drop its false self-citation; one misattribution fixed (`discrete_ordinates.rst` "V&V harness wiring" was never the degeneracy rule). Sphinx clean; grep=0.
- ~~**6 per-agent `AGENT.md` Tool-Freedom-Override sections**~~ **RESOLVED:** verified `.claude/rules/` auto-load into sub-agent contexts (inlined in the same project-instructions block as CLAUDE.md) → the sections were triple-redundant (rule + preloaded nexus skills + this) AND stale (their "default to Grep" framing is false since Grep/Glob were removed). Retired in all 6; agent-specific behavioral directives preserved under a slim `## Nexus`.
- ~~**Agent-index bloat**~~ **RESOLVED:** all 9 agents distilled — per-agent slim MEMORY.md index + behavioral `lessons.md` (indexes ~452 KB → ~64 KB; per-dispatch context tax eliminated). **git-authoritative throughout** (corrected ~7 SN campaigns mislabeled "in-flight"; only #236 open). qa/lessons deduped (L-034/039/040 + reorder); numerics L6-reorder + L8-pointer fixed.

**Also done this session (beyond the 3 findings):**
- **Main MEMORY.md** 27.2 KB → 4.5 KB, git-reconciled; 8 merged-campaign archaeology project memories retired.
- **4 meta-lessons graduated to `.claude/rules/`**: git-authoritative-merge-status, retirement-audit-completeness, mutation-test-hygiene, type-vs-property.
- **5 ceiling skills enriched**: coding-elegance (#19/#20 + Pattern-4/6), vv-principles (#12 + Mode-11 + Mode-10-doc-layer + offline-error), numerical-bug-signatures (Sig 8/9/10), cross-domain-frames (Smell #16 + A.2 + A.5), research (phantom-citation + catalogue≠method).
- **SN "Development history" changelog** seeded in `docs/theory/discrete_ordinates.rst` (17 git-verified milestones — the home for merged-campaign architectural history; #231 will distribute per-page slices).

**⚠ Topic-file retirement ATTEMPTED then REVERTED.** The finale `git rm`'d 333 merged-campaign closeouts — but they are NODES in the unified Nexus knowledge graph: **98+ are referenced from in-repo `plans/`+`docs/`** (incl. the SN theory page and the active pyright #226 keystone `issue_226_operator_generics_map.md`). The deletion was therefore NOT lossless — it violated the very **retirement-audit rule just graduated** (blast radius includes `docs/`+`plans/`, which `-W` won't warn on). Reverted (`git revert 0161129`): the closeouts are back as **COLD STORAGE** — referenced by the graph, NOT indexed by the slim MEMORY.md → **zero per-dispatch cost**. The per-dispatch win is the INDEX slim, NEVER file deletion. To delete them for real later: first clean their 98+ references (a separate careful sweep), or keep them as cold archive (recommended — they cost nothing loaded).

### (original flagged-findings text, now resolved, kept for the record)

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
