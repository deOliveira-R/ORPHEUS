# Instruction substrate & session-start context — findings and proposals

_Living decision document. Opened 2026-09-05 from the harness evaluation session; the raw chronological work log with every measurement is `scratch/_harness_eval/worklog_2026-09-05.md` (untracked — copy aside before any `git clean`). Drafts of every proposed surface: `scratch/_harness_eval/`._

**How to use this file.** Part I is what was MEASURED (`[M]` with the command or probe). Part II is what is PROPOSED, one entry per decision, each carrying the findings it rests on, its cost, its risk and a `Status`. Part III is the ledger of rulings. Part IV is the list of open questions to settle in discussion. Part V is the evidence appendix. Nothing in Part II has been implemented; nothing live was edited during the evaluation.

**Status vocabulary.** `PROPOSED` (not yet discussed) · `DISCUSSED` (talked through, no ruling) · `RULED <date>` (the user decided; the ruling is quoted) · `REFINED <date>` (ruled, then amended) · `REJECTED <date>` (with the reason, kept per plan-authoring §3) · `IMPLEMENTED @<hash>`.

**Goal.** The main agent and every sub-agent start with the same rules, lessons and V&V discipline they have today, organised by ROOT CAUSE rather than by order of discovery, enforced by tools where a check is mechanical, written for both Fable and Opus, at roughly a quarter of today's context cost.

---

# Part I — Findings

## F1. What loads at session start — `[M]` 2026-09-04

Token counts are `chars/3.6` (no tokenizer in the venv; the files are symbol-heavy, so true counts are likely higher). Commands: `wc -c`, a section-split script (worklog §1), one zero-tool haiku probe.

| block | who pays | ≈tokens | contents |
|---|---|---:|---|
| **A. always-on** | main agent AND every sub-agent dispatch | **64.6K** | plan-authoring.md 40.3K · coding-standards.md 9.6K · MEMORY.md 4.9K · process-discipline.md 3.3K · CLAUDE.md 2.7K · nexus-tools.md 2.5K · delegation.md 1.3K (+ vv-testing.md 0.9K, path-scoped `tests/**`) |
| **B. protocol batch** | main agent only (`session-start.txt` steps 1–6) | **≈125.5K** | vv-principles 42.5K + 2.6K `error_index.md` injected via `!cat` · lessons.md 40.9K · coding-elegance 21.6K · subagent-handoff 10.0K · `session_briefing()` ≈5.5K · development.rst 2.3K |
| **session start** | | **≈190K** | before the first user token, on top of the harness-fixed block |
| harness-fixed | main + every sub-agent | ≈12.9K + | 44 Nexus MCP tool schemas, eagerly loaded (`[M]` 43 visible to a haiku sub-agent) although the briefing's own `preload_hint` says they are deferred; plus 16 unrelated user-level skills (`moose-*`, `gitnexus-*`, `edgartools`, `graphify`) in every skill index |

Per-dispatch cost (A + AGENT.md + preloaded skills + agent-memory index): test-architect 168K · method-implementer 164K · numerics-investigator 159K · qa 156K · archivist 142K · elegance-enforcer 99K · literature-researcher 91K · cross-domain-attacker 86K · explorer 85K. **`[M]` empirical floor: a `general-purpose` haiku probe that answered ten yes/no questions with zero tool calls cost 113 006 tokens.**

Both Fable 5.1 (`[1m]` in user settings) and Opus 5 have a 1M window, so the cost is adherence and tokens per dispatch, not overflow.

## F2. Where the tokens sit, and how fast they arrived — `[M]`

| surface | composition |
|---|---|
| plan-authoring.md (1 022 lines) | 61 % is the surprise-log table (82 rows, 24.6K); **40 of 82 rows describe themselves as "no new clause" or a REPEAT**; §1–§10 = 15.0K with 21 `> [M]` worked cases; headings + first paragraphs = 1.2K; 119 `[M]` markers, 154 ⭐/⛔/⚠ glyphs |
| vv-principles (1 788 lines) | 36 anti-patterns = 24.6K (≈680 each; the NEVER sentence ≈40); test-design modes 11.6K; hierarchy/pillars/taxonomy/bit-identity ≈6K; 87 capitalised CRITICAL/NEVER/MUST tokens |
| lessons.md (2 475 lines) | 62 lessons; headlines = 1.7K (4 %); L1–L9 are 60–90 tokens each, L10+ average ≈700, L57 = 2 131 |
| coding-standards.md (451 lines) | 7.4K of 9.6K is one section ("Retire as you go") |
| coding-elegance (681 lines) | anti-patterns 4.4K; worked examples ≈1.2K; the 14-item checklist ≈1K |

Growth (`git log`, bytes at month end): plan-authoring **0 → 122K → 140K** (created 2026-08-06; 55 commits in August); lessons.md **72K → 149K** (Jul → Aug); vv-principles **47K → 138K → 153K** (Jul → Aug → Sep). Vendor guidance: "target under 200 lines per CLAUDE.md file. Longer files consume more context and reduce adherence"; SKILL.md "under 500 lines". Lesson cross-citation from rules/skills is small at the SYMBOL level (≈10 distinct lessons, after excluding V&V-level tokens `L0`–`L4`); the duplication is at the ROOT level (F5).

## F3. Who runs what — `[M]` 2026-09-05

Main agent: `claude-fable-5-1[1m]`. **7 of 9 project agents pin `model: opus`** (Opus 5); `explorer` and `literature-researcher` inherit Fable. Co-Authored-By trailers, last 400 commits: Fable 5 ×220, Opus 5 ×123, Fable 5.1 ×57. Roles (user): Fable and Opus do code development with user steering; **Haiku only for categorising fan-outs** over the whole tree; Sonnet not in use.

## F4. Harness mechanics, verified against the vendor docs and by probe

| mechanism | source | status |
|---|---|---|
| `.claude/rules/` is discovered **recursively** | memory doc | evidence must not live under `rules/` |
| path-scoped rules trigger **on Read** of a matching file | memory doc | usable for scoping |
| skills load `SKILL.md` only; `references/` on read; metadata (name + description) always in the index | skills docs | progressive disclosure available; `vv-principles/reference.md` (21.6 KB) already uses it |
| CLAUDE.md + rules + MEMORY.md are a **session-start snapshot** inherited by sub-agents; AGENT.md loads fresh per dispatch | memory 2026-06-22 + probe 2026-09-04 | substrate edits go live next session; agent edits next dispatch |
| MEMORY.md **is** inherited by sub-agents | `[M]` two probes (2026-06-22; 2026-09-05 quoted the first Active-work bullet verbatim) — **the doc says it is not** | keep it in the always-on budget |
| MEMORY.md auto-load cap: 200 lines / 25 KB | memory doc | today 63 lines / 17.8 KB |
| block HTML comments are stripped from CLAUDE.md at injection | memory doc | free maintainer notes; unverified for `rules/` |
| `InstructionsLoaded` hook logs which instruction files loaded and why; `/context` lists loaded memory files | memory doc | the acceptance instruments |
| **nested sub-agent dispatch works**: up to 3 layers below main by default (`CLAUDE_CODE_MAX_SUBAGENT_SPAWN_DEPTH`, `1` disables); `Agent` withheld at the depth limit or when an agent's `tools:` allowlist omits it; concurrency cap 20 (`CLAUDE_CODE_MAX_CONCURRENT_SUBAGENTS`) | sub-agents doc + `[M]` probe: depth-1 sonnet → depth-2 haiku → depth-3 haiku replied `PONG from depth 3`; depth-3 listed 65 tools, **no `Agent`** | the handoff skill's "Anthropic platform constraint: subagents cannot spawn" is DATED; the quoted sentence no longer exists on the cited page |
| sub-agents carry `SendMessage`; a named sibling can be messaged; a completed agent auto-resumes on message with full history | sub-agents doc + both probes listed `SendMessage` | the request/result relay is no longer the only path |
| per-agent `tools:` allowlists today | `[M]` frontmatter | `numerics-investigator`, `qa` include `Agent` (can dispatch now); the other 7 omit it |
| Nexus tool schemas eager in this harness | `[M]` 43 tools visible to a sub-agent | ≈12.9K per context, main + every dispatch |

## F5. Root-cause analyses — every clause is an instance of a small set of roots

Method: each mechanism/anti-pattern/pattern assigned by hand to the failure it instances (the full assignment tables are in the worklog §9–§10 and are meant to be disputed).

**plan-authoring.md — 96 mechanisms, 7 roots**

| root | one sentence | ≈mechanisms | log rows |
|---|---|---:|---|
| A. Description ≠ object | read the thing (precedent, pointer, blocker, memo, brief) before designing to its description | 10 | ~8 |
| B. Every claim carries its population and its instrument | denominator + predicate/filter validated per shape; numbers carry fixture, draw, protocol, precision | 23 | **~30 of 82** |
| B′. Authority does not transfer between neighbours | a marker, checkmark, positive control or rule vocabulary certifies only what its own command produced | 8 | ~7 |
| C. The campaign repeals its own facts | a plan is a snapshot; its own landings falsify tense; mark state where it is read | 11 | ~10 |
| D. A name is not an identity | pin a referent by index set, pairing, package, property | 6 | ~5 |
| E. An instrument must be able to fail | a gate/metric/canary/fix is evidence only if some realizable state changes its reading | 15 | ~11 |
| F. A step is a set on the dependence graph | every site that depends on the CONTRACT (most spell no symbol) + import edges both ways | 22 | ~17 |
| G. Outcome, not mechanism | goal in domain terms; means dated; done-when a predicate over the measured population | 7 | ~3 |

Diagnosis: the file is organised by ORDER OF DISCOVERY (§1…§10 are the dates their first case arrived), not by root. Root B is never named, so each new instance read as a new sharpening and was written as one — hence §2's 33 rows and the 40 self-described repeats. Five rows themselves conclude the check belongs at write-up time: that is a procedure, not a paragraph.

**vv-principles — 48 items, 5 roots**

| root | one sentence | ≈items |
|---|---|---:|
| V1. Independence is structural, not procedural | agreement is evidence only if no shared upstream (identity, integrand, input object, α-equivalent body) can produce it | 7 |
| V2. The instrument must be able to fail — (a) it runs, (b) the fault moves it on THIS fixture, (c) both proven by a positive control | the 6 test-design modes are (b); #17's seven sub-clauses are (c) | **21 of 48** |
| V3. Name the claim before choosing evidence | layer, level, quantity's definition, contract (type's invariant vs constructors' quality) | 8 |
| V4. Sample ≠ population; state the instrument | generated group, single draw, class inventory vs per-instance census, line vs window | 7 |
| V5. Bugs are term-local and sub-floor | six one-token substitutions; only per-term L0 sees them; assert the STRUCTURE the docstring names | 9 |

Diagnosis: V2 is never named; #1–#10 date from 2026-05, #11–#36 accreted through August one instance at a time — the same growth mechanism as plan-authoring §2.

**coding-elegance — 47 items, 5 roots**

| root | one sentence | ≈items | existing Nexus sweep |
|---|---|---:|---|
| C1. One quantity, one definition | value, convention, constant, primitive live at ONE site | 7 | `twin_paths` |
| C2. The type carries the invariant; discriminate once at the boundary | flags, strings, bare arrays, order-as-contract, deep re-checks are each a type never made | 12 | `discriminations`, `protocol_conformers` |
| C3. Code is notation for the math: name and dimension every quantity | | 10 | — |
| C4. Generality is decided by the domain's structure | primitive when algebra/symmetry/second instance says so; defer when only a name suggests it | 8 | `native_place` |
| C5. Debt compounds; clean now or write the trigger | | 4 | `dead_functions`, `dead_references` |

Diagnosis: the healthiest of the three — it already names one root explicitly ("a repeated conditional is a missing type; #3, #4, #7 are its corollaries") and its 14-item checklist is a procedure. Missing: the same move for the other four roots, and the root → tool → judgment structure that `nexus-elegance` carries in a separate file.

**Roots shared ACROSS files — the duplication a symbol grep cannot see**

| cross-file root | plan-authoring | vv-principles | elegance / standards / lessons |
|---|---|---|---|
| X1. The instrument must be able to fail | E | V2 | L44, L47, L54, L58; coding-standards demotion/promotion mirrors |
| X2. Sample ≠ population; state the instrument | B | V4 | L53, L55, L56, L59, L61, L62 |
| X3. Prose is not enforcement (docstring, marker, label, table header, plan claim ≠ what the code/gate asserts) | B′ | #14 #15 #20 #30 #33 #36 | elegance #20, #9; coding-standards tense + labelled-equation rules; L33 |
| X4. One definition per quantity | — | V1 | elegance C1; coding-standards retirement; L18, L21, L25 |

Each file grew by instances of X1–X3 arriving in whichever file the session was reading.

## F6. Vendor-guidance audit (pages fetched 2026-09-05: Prompting best practices; Prompting Claude Fable 5 / 5.1 / Opus 5; Skill authoring best practices; Claude Code memory & sub-agents)

**Both model pages endorse** (safe on shared surfaces): brief reasoned instructions over enumeration; the WHY beside the rule; examples in `<example>` tags; XML to separate content TYPES ("instructions, context, examples, variable inputs"); consistent terminology; scope-discipline text that both pages state in near-identical words; lead with the outcome.

**Where the pages diverge** and what it implies:
| topic | Fable 5 / 5.1 | Opus 5 | implication |
|---|---|---|---|
| verification nudges | "make self-verification explicit at intervals" | "remove — causes over-verification; avoid instructing re-checks it already performs" | generic nudges only on the Fable-only surface; DOMAIN V&V (mutation, positive controls, independence) is content and stays everywhere. `[M]` 3 such phrases in shared/agent surfaces, all domain-specific |
| review pre-filters | — | "only high-severity / be conservative" is followed literally → under-reports | never pre-filter in `qa`/`elegance-enforcer` briefs (`[M]` 0 today; matches adversarial-before-balance) |
| length | denser prose, fewer updates | longer responses AND longer files; calibrate explicitly | **`[M]` 0 of 9 AGENT.md files cap the return payload** — every Opus report lands unbounded in the main context |
| delegation | "use subagents frequently" | delegates readily; cap it | cap sentence on the two agents holding `Agent` |
| narration | ask for updates | tune down | main-agent only |
| prescriptiveness | "skills developed for prior models are often too prescriptive … can degrade output quality" | no statement; literal following is the risk | strip emphasis inflation on shared surfaces — safe on both |

**Practices the substrate is missing:** per-skill evaluations before documentation ("three scenarios, establish a baseline"); SKILL.md ≤ 500 lines; references one level deep with a TOC over 100 lines; no time-sensitive phrasing (`[M]` 5 dated-boundary phrases); enforcement in hooks rather than prose ("to block an action regardless of what Claude decides, use a PreToolUse hook"); one lesson per file with a one-line summary (Fable 5 memory guidance); consistent terminology.
**Checked clean:** no "show your reasoning" instructions (a Fable 5 refusal trigger); every skill/agent description < 1 024 chars, third person.

## F7. What the distilled drafts realised — `[M]` `scratch/_harness_eval/`

| draft | bytes | ≈tokens | preserved |
|---|---:|---:|---|
| plan-authoring.core.md | 25 567 | 7.1K | 96 mechanisms (appendix) — 4K target missed; merging would drop mechanisms |
| coding-standards.core.md | 15 976 | 4.4K | full retirement checklist + mirror + oracle exception |
| MEMORY.index.core.md | 8 247 | 2.3K | every entry; hooks ≤ 15 words |
| vv-principles.core.md | 37 515 | 10.4K | 36 + 6 + 6 items; ≈7.6K authored + 2.4K injected + split table |
| coding-elegance.core.md | 32 586 | 9.1K | 7 patterns, 20 anti-patterns, 14-item checklist verbatim |
| subagent-handoff.core.md | 13 464 | 3.7K | table + triggers + 3 block contracts; main half 2.1K / sub half 0.9K |
| lessons.index.md | 12 540 | 3.5K | 63 lines; 18 fully absorbed, 22 headline-hides-check, 2 archaeology |
| session-start.txt.proposed | 1 694 | 0.5K | same gates, budgeted batch |
| plan-checklist.md | 4 774 | 1.3K | floor (11) + ceiling |

Realised: **tier 1 ≈ 23.6K** (vs 64.6K; process-discipline, nexus-tools, delegation, CLAUDE.md untrimmed), **tier 2 ≈ 32.2K** (vs 125.5K), **session start ≈ 56K vs ≈190K**; no-op dispatch floor ≈ 70K vs 113K.

## F8. The docs → eager-index pattern already exists for one surface — `[M]` 2026-09-05

`tools/verification/generate_error_index.py` reads `docs/theory/verification/error_catalog.rst` (488 018 B, 82 `.. error-entry::` entries, 293 catching tests) and regenerates `.claude/skills/vv-principles/error_index.md` (9.3 KB), whose header reads *"Bodies live once, in the corpus. This index is derived from the graph; editing it by hand is a no-op."* It has a `--check` mode for CI and runs on every Sphinx build. The catalogue entries are Nexus nodes (`vv:error:ERR-NNN`), so `errors()`, `dead_references` and `staleness` see them. This is exactly the SKILL.md-eager / reference-lazy split with the DOCS as the reference layer — proven on the largest single knowledge surface in the repo.

`docs/development.rst` (8.4 KB, git workflow) is a single page in the root toctree; there is no `docs/development/` section yet. Agents already read and write RST daily.

## F9. Markdown in the docs, and typed blocks — `[M]` 2026-09-08

- **MyST is available and compatible:** `pip install --dry-run myst-parser` on this venv resolves to **myst-parser 5.1.0 + markdown-it-py 4.2.0** against Sphinx 9.1.0, nothing else touched. It is not installed today; `docs/` carries no `.md` file.
- **MyST reaches the whole Sphinx directive/role set from Markdown**: fenced ```` ```{directive} ```` / `:::{directive}` blocks with `:option:` lines; `{role}` inline; `(label)=` targets; `{eval-rst}` for raw RST. So a Markdown page can carry `{error-entry}`, `{verifies}`, `{math}` — the same typed blocks as RST.
- **An RST page can include a Markdown file parsed by MyST**: docutils `include` has a `:parser:` option ("Parse the included content with the specified parser", docutils ≥ 0.17, provisional) — the MyST-documented idiom is `.. include:: ../README.md` + `:parser: myst_parser.sphinx_`. A rule or skill core that must stay `.md` for the harness can therefore be RENDERED in the docs from its single source.
- **Nexus is format-agnostic:** it registers 8 typed directives (`error-entry`, `verifies`, `implements`, `no-implementation`, `discretizes`, `derives-from`, `approximates`, `nexus-graph`) at `doctree-read` and extracts from `env.get_doctree()` — sections, labels, `pending_xref`, math. A MyST page produces the same doctree, so it is indexed identically.
- **The corpus already has three ad-hoc "typed block" conventions**: `.. admonition::` with `:class:` (100 uses: "Key Facts", "Development history …"), `sphinx_design` `.. dropdown::` (12: the `nexus-meta` machine headers), and **470 `.. (vv-status rationale)` RST COMMENTS** — a tagging convention that is invisible to the build and to Nexus, which is exactly why the 2026-08-19 surprise ("2 of 94 equations carry one") had to be answered by grep.
- **The docs already cite archived Markdown plans as authoritative** (`loss_representation.rst:3792` *"`.claude/plans/archive/sn_sweep_strategy.md` — the authoritative …"*, `diffusion_1d.rst:1447`, …) by PATH in a literal — unrenderable, unvalidated, unindexed. Durable Markdown is already in the corpus's dependency graph; it just is not in the corpus.
- `.. only::` is a BUILD-tag mechanism ("control only content of document" via `-t tag`) — conditional inclusion, not semantic typing. `.. container::` / `.. class::` set CSS classes only: no options, no validation, no node type.

## F10. `omitClaudeMd: true` drops CLAUDE.md, EVERY rule file, AND the auto-memory index — `[M]` 2026-09-20, Claude Code 2.1.278

The sub-agents doc documents a per-agent frontmatter field **`omitClaudeMd: true`** (*"launch this subagent without the user, project, and local CLAUDE.md files; managed policy files still load … Requires Claude Code v2.1.271 or later"*). On 2026-09-19 this machine ran 2.1.261 and the field was inert; the user upgraded to **2.1.278** and the probe ran 2026-09-20.

**Method.** Two throwaway agents, identical except for the flag (`model: haiku`, one-line system prompt, no tools used), dispatched by a headless `claude -p --model sonnet` session (the running session's agent registry is fixed at start, so new agent directories are invisible until a fresh session) with the same 10-question zero-tool prompt (`scratch/_harness_eval/probe_questions.md`); usage read from the sub-agent transcripts under `~/.claude/projects/…/c03a470c…/subagents/` (0 `tool_use` blocks in both). The throwaway agents were removed afterwards; `git status` on `.claude/agents/` is clean.

| context item | keep (control) | **omit** |
|---|---|---|
| CLAUDE.md heading | YES | **NO** |
| `rules/plan-authoring.md` (surprise-log heading) | YES | **NO** |
| `rules/coding-standards.md` ("Retire as you go") | YES | **NO** |
| `rules/nexus-tools.md` ("ugrep 7.5.0") | YES | **NO** |
| auto-memory `MEMORY.md` ("Four index disciplines" + the 60 chars after it; first Active-work bullet) | YES, quoted verbatim | **NO** |
| Co-Authored-By trailer instruction | YES | YES (harness-injected, not from CLAUDE.md) |
| SessionStart hook text (`ENVIRONMENT HEALTH`) | NO | NO (hooks do not reach sub-agents) |
| Nexus MCP tool schemas | ~46 | ~45 (unchanged — the harness-fixed block) |
| **first-turn input tokens** (`cache_creation + cache_read`) | **107 476** | **36 510** |

So the switch removes **≈71K tokens** from a dispatch — CLAUDE.md, the whole `.claude/rules/` set (loaded "with the same priority as `.claude/CLAUDE.md`", and dropped with it) and MEMORY.md together — leaving the harness-fixed block (tool schemas ≈13K, the skill index, the agent roster, the system prompt) at ≈36.5K. That is a larger per-dispatch saving than the entire A1 distillation delivers (64.6K → 23.6K), and it is available NOW for every agent that can take its rules from the brief. ⚠ Two consequences: (a) an `omitClaudeMd` agent sees NO project rule at all — not the ugrep hazard, not the canonical pytest spelling, not the local-folder-first rule — so its brief template becomes its only rule surface (K8); (b) it also loses MEMORY.md, which VII.3 showed IS inherited otherwise, so a Support agent that needed a campaign pointer must be handed it. The previously measured 113K floor for a `general-purpose` haiku probe is consistent: 107.5K here on a smaller system prompt.

⚠ Side finding from the same run: the harness now warns *"Permission allow rule (.claude/settings.json): `Write(.claude/agents/**)` is not matched by file permission checks — only `Edit(path)` rules are. Use `Edit(.claude/agents/**)` instead"* — a one-line settings fix (Group C housekeeping).

---

# Part II — Proposals

Each entry: **What** · **Why** (findings it rests on) · **Where** · **Cost / risk** · **Status** · **Open**. IDs are stable; cite them in rulings.

## Group A — Structure of the substrate (load-time tiers)

**A1. Tier every instruction surface by load time.** Always-on core (main + every dispatch) ≤ ~24K; main-agent protocol batch ≤ ~32K; everything else on demand behind a pointer. — Why: F1, F2, F7. — Where: all of `.claude/`. — Cost: one-time restructure (P1–P3 in G8); every rule keeps imperative + check + tell + dated pointer, so nothing is lost, only moved. Risk: a clause whose persuasion lived in its story loses force as a one-liner; mitigated by the `tl:` field and detected by G5. — Status: PROPOSED.

**A2. An evidence directory for founding cases.** `.claude/evidence/<rule>.md` holds every `> [M]` blockquote and the surprise-log table verbatim, under anchors the cores point at; `.claude/evidence/lessons.md` holds lesson bodies under `## L<n>`. — Why: F4 (`rules/` is recursive, so evidence cannot live there); F2. — Cost: a mechanical move script + an anchor checker (G3). — Status: SUPERSEDED by I1/I4 if Q17 rules for the docs (2026-09-05). — Open: Q1 (location vs HTML-comment trick).

**A3. Skills become core + `references/`.** vv-principles, coding-elegance, subagent-handoff first; phase 2: numerical-bug-signatures (10.1K, preloaded by 4 agents) and algebra-of-record (9.0K, by 3). — Why: F1 (skills are 60 % of the batch and 55–95K of every heavy dispatch), F4 (progressive disclosure), F6 (≤ 500 lines). — Cost: the drafts exist (F7); references are the moved text. Risk: none in-session — skill edits are live at next invocation, so this phase can be smoke-tested immediately. — Status: PROPOSED.

**A4. lessons.md becomes a one-line index** (imperative · check · twin); bodies to `evidence/lessons.md`; the 18 lessons fully absorbed by a rule/skill clause keep only their `twin:` pointer; L6 (a physics fact) moves to its theory page; L21 is archived with its superseded banner. — Why: F2 (headlines = 4 %), F5 (X1–X4 duplication), F6 (one lesson per file with a one-line summary is the Fable 5 memory guidance verbatim). — Status: PROPOSED; the bodies' home moves to `docs/development/lessons.rst` under I1. — Open: Q3 (which twin survives).

**A5. MEMORY.md re-indexed by its own disciplines**, plus a fifth: a hook is ≤ 15 words. Draft 8.2 KB vs 17.8 KB. — Why: F1, F4 (inherited by sub-agents; 25 KB cap). — Cost: outside git — back up to `scratch/` first; live next session. — Status: PROPOSED.

**A6. session-start.txt: same gates, budgeted batch.** Cores + lessons index + briefing; `development.rst` on demand; the batch states its token budget so drift is visible. — Why: F1. — Status: PROPOSED.

**A7. AGENT.md bodies slimmed.** `archivist` (30.8 KB) and `elegance-enforcer` (34.9 KB) carry 8–9 KB "institutional knowledge" sections that belong in their agent-memory topic files. — Why: F1. — Status: PROPOSED.

**A8. Harness-side items (user settings, not the repo).** (a) Nexus tool schemas: 12.9K eager per context; the briefing believes they are deferred — worth raising with the nexus repo / harness setting. (b) The 16 unrelated user-level skills could move to a per-project plugin. — Why: F1, F4. — Status: PROPOSED.

## Group B — Reorganise content by ROOT, not by order of discovery

**B1. plan-authoring.md organised by roots A–G.** Each root: one sentence, its checks, its tells, pointers to evidence. Keep the existing § numbers as stable IDs inside the root layer (they are cited from 82 rows, memories and issues). — Why: F5 (§2 = 33 rows because root B was never named). — Risk: renumbering breaks citations — hence keep §-IDs. — Status: PROPOSED. — Open: Q5.

**B2. vv-principles organised by V1–V5.** The six test-design modes and #17's seven sub-clauses fold under V2 (a)/(b)/(c); item numbers kept as IDs. — Why: F5 (V2 = 21 of 48, never named). — Status: PROPOSED.

**B3. coding-elegance organised by C1–C5 as root → run the Nexus sweep → judge the residue**; consider merging with `nexus-elegance` (4.1K), which is that mapping today in a second file. — Why: F5. — Status: PROPOSED. — Open: Q7.

**B4. A shared instrument doctrine** (X1 the instrument must be able to fail; X2 sample ≠ population, state the instrument; X3 prose is not enforcement; X4 one definition per quantity): one page, cited by plan-authoring, vv-principles, coding-elegance, coding-standards; each file keeps only its domain-specific instances. — Why: F5 (the same four roots spelled in four vocabularies). — Status: PROPOSED. — Open: Q8 (rule vs skill).

**B5. Checklists as the floor.** Plan floor (11 items, `plan-checklist.md`) + ceiling; verification-plan floor (7 items, worklog §10d); elegance checklist stays (14). Each item names its root. — Why: user request; F5. — Status: PROPOSED.

**B6. The surprise log keeps appending, tagged by root.** A row records which root (A–G / V / C / X) it instances; a root that recurs is the trigger for a TOOL (Group C), not a new paragraph. The file's own target (surprises → 0) becomes measurable per root. — Why: F2, F5. — Status: PROPOSED.

## Group C — Tools and hooks where the check is mechanical

**C1. `tools/plan/claim_lint.py`** (root B / X2): over a plan, issue body or commit-message diff, flag every universal (every/all/none/each/only/no) without `of N` in the sentence, every bare number without `[M]` + a command within two lines, every ratio with one population. — Why: F5 (~30 of 82 rows), F6 (hooks over prose). — Status: PROPOSED. — Open: Q9 (also as a hook?).

**C2. `tools/plan/blast_set.py <symbol|class|alias>`** (root F): AST census by CONTRACT — literal calls, `__subclasses__` at runtime, registry loops, `getattr("…")` strings, alias expansion + comparison partners, einsum subscript letters, import edge both ways with relative imports resolved and intra-file namers; prints its own predicate limits. The 13 §6b spellings become its test suite. — Status: PROPOSED.

**C3. `tools/plan/date_collision.py`** (root C): for each `[M] <date>` in a plan, list campaign landings on/after it and mark tree-half claims RE-RUN. — Status: PROPOSED.

**C4. Activation-counter fixture** (root E / X1): the spy used at P4.9a made reusable — "does this canary execute the carved path at all?" — Status: PROPOSED.

**C5. Hooks for enforcement-shaped rules.** PreToolUse: deny `git add -A`; guard `git checkout/restore/stash` on a path with uncommitted edits. PostToolUse: after a rename/delete, remind `dead_references`. `InstructionsLoaded`: log what loaded (G4). — Why: F6 ("use a PreToolUse hook"); each costs zero context. — Status: PROPOSED.

**C6. Elegance roots run their Nexus sweeps first** (`twin_paths`, `discriminations`, `protocol_conformers`, `native_place`, `dead_functions`, `dead_references`) — already-built tools, currently referenced from a separate skill. — Status: PROPOSED (folds into B3).

## Group D — Model-aware surfaces (Fable main agent, Opus sub-agents, Haiku fan-outs)

**D1. Shared surfaces carry only what both model pages endorse** (F6 first paragraph). — Status: PROPOSED.

**D2. Fable-only nudges live in `session-start.txt`** (read by the main agent only): finish the whole task, progress updates, explicit verification cadence for long runs. — Why: F3, F6 (verification-nudge conflict). — Status: PROPOSED.

**D3. Opus-pinned AGENT.md files get:** a return-payload cap ("report under N words; the file carries the detail"); no generic re-check nudges; no severity pre-filters in review briefs; on `qa` and `numerics-investigator` the Opus-page delegation sentence. — Why: F6 (`[M]` 0 of 9 carry a cap; Opus 5 runs long; delegates readily). — Status: PROPOSED.

**D4. Haiku fan-out brief template**: low-freedom, schema-fixed instructions (skill doc: "does the Skill provide enough guidance?"); for headless fan-outs (`claude -p`) use `--setting-sources` without `project` to skip the always-on block entirely. — Why: F3, F1 (the 65K block is pure cost to a categoriser). — Status: PROPOSED.

**D5. Effort per agent** where the harness exposes it (census/explorer tasks at lower effort; Opus 5: "use low and medium liberally wherever quality holds"). — Status: PROPOSED (needs a harness check).

## Group E — Sub-agent handoff protocol update

**E1. Retire the "Anthropic platform constraint" claim** _(stands as H4's first step)_ in `subagent-handoff-protocol/SKILL.md:17-24` and `CLAUDE.md:102`; replace with the real mechanism: per-agent `tools:` allowlist decides; depth cap (default 3); concurrency cap (default 20). — Why: F4 (doc + `PONG from depth 3`). — Status: PROPOSED.

**E2. Decide per agent whether it SHOULD spawn.** Recommendation: keep `Agent` only on `numerics-investigator` and `qa` (probe fan-outs, verifier patterns); set `CLAUDE_CODE_MAX_SUBAGENT_SPAWN_DEPTH=2` in `.claude/settings.json` `env` so a grandchild cannot spawn; leave concurrency at 20. — Why: F4, F1 (each nested dispatch re-pays the always-on block one level down), F6 (Opus 5 delegates readily). — Status: SUPERSEDED by Group H (2026-09-05), kept per plan-authoring §3. — Open: Q6.

**E3. `SendMessage` becomes the coordination primitive**; the DISPATCH_REQUEST / DISPATCH_RESULT blocks stay as the fallback contract for agents without `Agent`, moved to `references/`; the "capture agent_id" section is superseded by names. — Why: F4. — Status: SUPERSEDED by Group H (2026-09-05), kept per plan-authoring §3.

**E4. Split the skill's two halves**: main-agent half (agent table + triggers, 2.1K) preloaded by nobody but read at session start; sub-agent half (block formats, 0.9K) preloaded by the 8 agents that list the skill today. — Why: F7 (the split is clean). — Status: SUPERSEDED by Group H (2026-09-05), kept per plan-authoring §3.

## Group F — Writing conventions

**F1. XML tags for CONTENT TYPES, not sections**: `<evidence clause="§1" date="…">`, `<example>`, `<contract>` (the dispatch blocks), `<summary>` (compaction). Fixed vocabulary across files; Markdown keeps structure; no tags in frontmatter `name`/`description`; ⚠ `<name>`, `<path>`, `<hash>` are already PLACEHOLDERS (`[M]` 35/19/6 hits) — no collision. — Why: F6 (the doc's stated case: "instructions, context, examples, variable inputs" mixed). — Status: PROPOSED.

**F2. Emphasis reduction.** Rule + why + check + tell; ⛔ reserved for a live refutation; ⭐/⚠ retired from the cores. — Why: F2 (154 glyphs, 87 caps), F6. — Status: PROPOSED. — Open: Q10.

**F3. Consistent terminology**: one term per concept (gate/pin/witness; floor/ceiling/core) chosen in the cores. — Why: F6. — Status: PROPOSED.

**F4. Time-sensitive phrasing** → "old patterns" sections (`[M]` 5 hits). — Status: PROPOSED.

**F5. HTML comments for maintainer notes** in CLAUDE.md (doc-verified) — test for `rules/` with the `InstructionsLoaded` hook before relying on it. — Status: PROPOSED.

## Group G — Acceptance, validation, sequencing

**G1. `/context` in a FRESH session** before/after: memory-files block ≤ 25K (today ≈ 65K), batch ≤ 35K (today ≈ 125K). Rules/CLAUDE.md/MEMORY.md edits are a session-start snapshot, so nothing can be verified in the editing session.
**G2. No-op sub-agent probe** ≤ 75K tokens (today `[M]` 113 006).
**G3. Anchor check**: every `→` pointer in a core resolves to a heading in its evidence file (0 dead); mechanism-count appendix = source census.
**G4. `InstructionsLoaded` hook log** confirming exactly which files loaded, for main and for a sub-agent.
**G5. Surprises per campaign, tagged by root**, on the first campaign step executed under the new substrate (CS4c step 5).
**G6. A/B on BOTH models** (skill doc: test with every model you use): each core once as a Fable main-session load, once as an Opus sub-agent preload; optional `qa` review of the step-5 diff with full vs core vv-principles — n=1, a smoke test.
**G7. Per-skill evaluations** (three scenarios + baseline) — the vendor practice we do not have; the tool-routing ablation plan is the seed.
**G8. Sequencing (each phase one commit, gated):** P1 evidence dir + rules split + anchor check → P2 skills split (smoke-testable in-session) → P3 lessons index, MEMORY.md, session-start.txt, CLAUDE.md trims → P4 acceptance in a fresh session → P5 agent preloads and AGENT.md bodies → P6 tools and hooks (C1–C5) → P7 root reorganisation (B1–B4), which can proceed independently of P1–P3 once the ordering question (Q5, Q11) is ruled.
— Status of G1–G8: PROPOSED.

## Group H — Workflows and agent roles (RESHAPES Group E: the handoff protocol was built for a world without nested dispatch)

**Premise (user, 2026-09-05).** Inter-agent communication is now a tool grant, not a protocol; the remedy is to encode the WORKFLOWS we want and the roles inside them — routing guidance, explicitly non-exhaustive, with agents dispatching whenever the reason is worthwhile.

**H1. Role taxonomy.** Three roles, and an agent may hold two.
| role | mandate | `Agent` tool | agents |
|---|---|---|---|
| **Orchestrator** | owns phase transitions, the review stage, user rulings, issues, commits | yes (depth 0) | the main agent |
| **Key** | owns one phase of a workflow; may call Support agents | yes | test-architect · method-implementer · numerics-investigator · archivist · qa |
| **Support** | stateless helper: answers a question and returns; never spawns | no (allowlist) — and withheld by the depth cap anyway | explorer (available to ANY agent) · literature-researcher (papers, local folder first, OCR) · cross-domain-attacker (adversarial detection of native structure) · haiku categorisers (`general-purpose`, schema-fixed briefs) |
| **Review** | judges work it did not author; dispatched by the Orchestrator only | qa: yes (for Support calls); elegance-enforcer: no | qa · elegance-enforcer · cross-domain-attacker (in its post-first-pass role) |
— Why: F4 (nesting works; `tools:` decides), F1 (each nested dispatch re-pays the always-on block), F6 (Opus 5 delegates readily; "do not use subagents to verify your own work"). — Status: DISCUSSED (the user named the roles; the table is the proposed encoding). — Open: Q14.

**H2. Hard invariants (the only fences; everything else is guidance).**
1. **Authors never dispatch their own reviewers.** qa and elegance-enforcer are dispatched by the Orchestrator on the artefact, never by the implementer — structural independence (vv-principles V1) and the Opus 5 page in one rule.
2. **Depth ≤ 2**: Orchestrator → Key → Support. Set `CLAUDE_CODE_MAX_SUBAGENT_SPAWN_DEPTH=2`; Support allowlists omit `Agent`.
3. **Continuity by name.** The test-architect that wrote the verification spec is RESUMED (`SendMessage`) at review time to confirm its gates landed; the implementer is resumed with the review findings. No re-briefing, no lost context.
4. **Every brief carries** the workflow ID, the phase, the previous phase's artefact paths, and the return contract: a word cap (Opus runs long), files carry the detail, and a `NEEDS:` block for anything the agent could not obtain (the minimal descendant of DISPATCH_REQUEST).
— Status: PROPOSED. — Open: Q15.

**H3. The workflows** (phase → key agent → supports → gate → next). Guidance, not a fence.
| ID | workflow | phases |
|---|---|---|
| **W1** | Build a capability (feature / method) | P0 context: explorer (+ literature-researcher when the formulation is published) → P1 verification design: **test-architect** (spec + gates naming their first red) → P2 build: **method-implementer**, or the main agent for surgical carves (delegation.md exception) — supports explorer, literature-researcher; cross-domain-attacker after the first pass → P3 review, PARALLEL, dispatched by the Orchestrator: **qa** (term-level correctness, coverage, mutation) + **elegance-enforcer** (structure; three-leg VIOLATION verdicts) → gate: any red ⟹ resume the implementer with the findings → P4 documentation: **archivist** (theory page, derivations, changelog; gates: `sphinx -W`, `dead_references` 0) → P5 close-out: Orchestrator (issues, retirement audit, commit) |
| **W2** | Wrong answer / bug | P1 **numerics-investigator** (probe cascade; supports explorer, literature-researcher; may call test-architect for the permanent test) → P2 fix (implementer or main) → P3 **qa** (fix + mutation) → P4 **archivist** (ERR-NNN entry + theory page) → close-out |
| **W3** | Surgical carve / refactor (user-steered) | main agent writes; explorer for the blast set (later C2); **test-architect** for gates and re-baselines; review as W1-P3; **archivist** for the changelog |
| **W4** | Documentation campaign | **archivist** key; qa for claim verification; explorer; gates `-W`, `dead_references`, `staleness` |
| **W5** | Design / adversarial review | **cross-domain-attacker** + **elegance-enforcer** on a first-pass design; output feeds W1-P2 |
| **W6** | Tree-wide census / fan-out | haiku `general-purpose` categorisers with a fixed output schema (D4); the Orchestrator aggregates; no nesting |
| **W7** | Literature acquisition / OCR | **literature-researcher** (support): `scratch/literature/` first, OCR sidecars, "not in local folder" is a question to the user, never a pivot |
— Status: PROPOSED (W1 is the user's stated known workflow; W2–W7 are inferred from the agents' mandates).

**H4. Encoding.** (a) `.claude/rules/workflows.md`, always-on, ≤ 2K tokens: the roles table, the four invariants, the seven workflow one-liners — this replaces the agent table + proactive triggers half of the handoff skill (2.1K). (b) Each AGENT.md gains a ten-line "my phase, my supports, my return contract" block (live at next dispatch). (c) The full phase descriptions, brief templates and return contract live in `docs/development/workflows.rst` (Group I), read on demand. (d) `subagent-handoff-protocol` is RETIRED; its block formats survive only as the `NEEDS:` return contract. — Status: PROPOSED. — Open: Q16.

**H5. Frontmatter changes implied.** Add `Agent` to test-architect, method-implementer, archivist; keep it on numerics-investigator and qa; keep it OFF explorer, literature-researcher, cross-domain-attacker, elegance-enforcer. Add the Opus delegation sentence to every Key agent. — Status: PROPOSED.

## Group I — `docs/development/` as the on-demand knowledge layer (RESHAPES A2, A4, part of A3)

**Premise (user, 2026-09-05).** Open a "development" section in the docs to hold what agents should be able to consult without pre-loading it — the SKILL.md-eager / reference-lazy split, with the corpus as the reference layer. F8 shows this already works for the error catalogue.

**I1. The section.** `docs/development/index.rst` (today's `development.rst` becomes `git_workflow.rst` inside it) with pages: `workflows.rst` (H3–H4 in full) · `plan_authoring.rst` (every founding case + the surprise log, one labelled section per row: `.. _surprise-2026-08-14-precedent:`) · `lessons.rst` (one labelled section per lesson) · `vv_anti_patterns.rst` + `test_design_modes.rst` (mechanisms + founding cases) · `elegance_patterns.rst` (worked examples, code contrasts) · `retirement_audit.rst` (coding-standards evidence) · `harness.rst` (what loads when, the measured budgets, how to add a rule/skill/agent — this evaluation's durable residue). — Why: F8, Cardinal Rule 3 (Sphinx IS the brain — extended from theory to process knowledge). — Status: DISCUSSED. — Open: Q17.

**I2. The eager layer points into the docs, and the INDEX files are GENERATED from them** _(= #308 channel 2; see Part VI)_. Generalise `generate_error_index.py` → `tools/docs/generate_indexes.py`: the lessons one-line index from `lessons.rst` (title + `:check:` + `:twin:` fields), the vv anti-pattern one-liners from `vv_anti_patterns.rst` (`:never:` / `:instead:` / `:check:` / `:tell:`), the plan-authoring core's pointer table from `plan_authoring.rst`. Bodies live once; `--check` in CI; drift impossible (elegance C1). — Why: F8, X4. — Status: PROPOSED. — Open: Q18 (directives vs plain sections).

**I3. Gates that come for free.** `sphinx -W` validates every `:ref:` (G3's anchor check becomes the docs build); Nexus `dead_references` catches a retired symbol still cited in a founding case (today the surprise log names `AngularAverageOperator`, `SweepCoefficientCache` … and nothing checks them); `staleness` flags evidence whose cited code moved; the nexus PostToolUse hook already fires on every edit. — Status: PROPOSED.

**I4. What stays in `.claude/`.** Only harness-facing loaders: CLAUDE.md, rule cores, skill cores, agents, hooks, memory, and plans (transient; a plan's compaction record feeds `docs/development/` at close-out, not before). `.claude/evidence/` (A2) is NOT created if Q17 rules for the docs. — Status: PROPOSED. — Open: Q19.

**I5. Migration hazards, named up front.** Markdown → RST conversion of ~150 KB of surprise log + blockquotes: the 2026-08-19 backtick-mangling case applies — use a CONSERVATION check (count of inline literals in ≡ out), not a warning count; per-row sections rather than a 143-char/line table so `:ref:` labels are per row; the corpus grows ~600 KB (the error catalogue is already 488 KB, so build time is not the constraint). — Status: PROPOSED.

**I6. Markdown vs RST — the criterion (RULED premise 2026-09-08: plans are transient and stay in `.claude/plans/`; the docs are durable; no `.claude/evidence/`).**
A durable file is RST unless a SECOND CONSUMER requires Markdown. Test: *"who else reads this file, and in what format?"*
| durable content | format | why |
|---|---|---|
| theory, derivations, equations, catalogue entries, anything owning `:eq:` / `:cite:` / typed directives, anything the archivist authors | **RST** | authored for the docs; the Sphinx domain machinery is the point |
| rule cores, skill cores, agent definitions — harness-loaded, MUST be `.md` | ~~MD, single-sourced via include into the docs~~ **[REFUTED 2026-09-19 by the Q20 ruling: the flow is docs → `.claude`, never the reverse]** → SOURCE in `docs/development/` (MyST), `.claude/` copy GENERATED and committed (K7) | the docs are the harness-independent brain; `.claude/` is one harness's view of it |
| campaign plans at close-out (compaction record, rulings, memos), archived plans the theory pages already cite | **MD, moved verbatim** to `docs/development/archive/` | authored in MD during work; conversion is where content is lost (the 2026-08-19 backtick case); the docs then VALIDATE and INDEX what they cite today by bare path |
| founding cases / surprise log / lesson bodies / anti-pattern mechanisms (the evidence layer) | ~~RST, written fresh~~ → **MD (MyST) with `:::{evidence}` fences** — amended 2026-09-08 by #422 (Part VI.3) | dense nested inline markup that RST cannot express; a move, not a conversion |
| README, issue templates, GitHub-facing text | MD | GitHub is the consumer |
Typing does not depend on format: an MD page carries the same blocks as `:::{evidence}` fences. — Status: PROPOSED (criterion); the premise is RULED. — Open: Q20.

**I7. Typed blocks: Sphinx directives, not XML tags** _(= #308 channel 3 with typed fields; validation depends on nexus#90's fix, Part VI)_. What the vendor doc wants from XML tags (unambiguous separation of instruction, example, evidence, context) is delivered in the docs by DIRECTIVES, which do three things a raw `<evidence>` tag cannot: the build validates them (`-W`: unknown directive or missing required option is an error), Nexus turns them into graph nodes, and they render for humans — while in SOURCE, which is what the model reads, `.. evidence::` / `:::{evidence}` is exactly as explicit as a tag. Proposed set, registered by a small `docs/_ext/devblocks.py` (or contributed to `sphinxcontrib-nexus`, which already owns the pattern):
| directive | required fields | graph node / edge | replaces |
|---|---|---|---|
| `.. rule:: <id>` | `:root:` (A–G / V1–V5 / C1–C5 / X1–X4), `:check:`, `:tell:` | `dev:rule:<id>` | the clause bullets of the cores |
| `.. evidence:: <id>` | `:date:`, `:rule:` | `dev:evidence:<id>` → `instances` → rule | every `> [M]` blockquote and surprise-log row |
| `.. example::` | `:kind:` positive / negative / counter | child of its rule | worked examples in skills |
| `.. lesson:: L<n>` | `:check:`, `:twin:` | `dev:lesson:L<n>` | lessons.md bodies |
| `.. anti-pattern:: <id>` | `:never:`, `:instead:`, `:check:`, `:tell:` | `dev:antipattern:<id>` | the vv and elegance catalogues |
| `.. rationale::` | `:status:` | `dev:rationale` → the equation it follows | the 470 `.. (vv-status rationale)` comments — the question "which equations have one?" becomes a Nexus query |
Rendering: each is an admonition with a class (the `Key Facts` idiom the corpus already uses), so nothing changes for a human reader. The eager index files (`I2`) are generated FROM these nodes, the way `error_index.md` is generated from `error-entry` nodes today. Rejected alternatives, with reasons: raw XML tags in RST (inert text: no validation, no node, no rendering); `.. only::` (build-conditional, not semantic); `.. container::` / `.. class::` (CSS only, no fields); plain `.. admonition:: Evidence` (works today with zero extension code and is the interim fallback, but carries no fields and mints no node). — Status: PROPOSED. — Open: Q18 (now concrete), Q21.

**I8. Inline marker roles (optional, later).** `:M:` (measured: `` :M:`2026-09-04, wc -c` `` renders `[M]` with its command), `:R:` (reasoned), `:refuted:`. Makes plan-authoring's marker vocabulary machine-checkable in the docs — "a claim with no command" becomes a query, the C1 lint's docs-side half. — Status: PROPOSED (deferred until I7 lands).

---

# Part III — Decision ledger

| ID | status | ruling (quoted) | date |
|---|---|---|---|
| A1–A8 | PROPOSED | | |
| B1–B6 | PROPOSED | | |
| C1–C6 | PROPOSED | | |
| D1–D5 | PROPOSED | | |
| E1–E4 | PROPOSED | | |
| F1–F5 | PROPOSED | | |
| G1–G8 | PROPOSED | | |
| H1–H5 | H1 DISCUSSED, rest PROPOSED | user named the roles and W1 (2026-09-05) | 2026-09-05 |
| I1–I5 | I1 DISCUSSED, rest PROPOSED | user proposed the docs section (2026-09-05) | 2026-09-05 |
| E2–E4, A2 | SUPERSEDED | by H / I | 2026-09-05 |
| Q19 → I4 | RULED | "Plans should stay in .claude/plans/ because they are transient files. We execute them, triage, then archive (which is a staging spot for exclusion). Docs are durable files." | 2026-09-08 |
| Q17 → I1 | RULED (home) | "An evidence folder in .claude/ goes exactly against this point" — durable evidence lives in the docs; the MD-vs-RST question is I6 | 2026-09-08 |
| I6–I8 | PROPOSED | user asked to check MyST and Sphinx tags first (2026-09-08); checked: F9; I6 amended by #422 (Part VI.3) | 2026-09-08 |
| J1–J5 | PROPOSED | issue sweep 2026-09-08 (Part VI) | 2026-09-08 |
| Q1,2,3,5,6,7,9,11,14,15,16,18,20,21 | RULED | Part VII.1, verbatim where the words carry the criterion | 2026-09-19 |
| Q12 | SETTLED | inherited (0 tool calls, verbatim continuation) | 2026-09-19 |
| H1, H2.1, H5 | REFINED | every key agent spawns; elegance-enforcer is key; parent review independent | 2026-09-19 |
| I6 (rule-core row) | REFUTED → K7 | flow is docs → .claude only | 2026-09-19 |
| K1–K7, J6 | PROPOSED | this round | 2026-09-19 |
| Q4, Q8, Q13, Q10, Q22, Q23, Q24 | RULED | round 2: recall test; rule+skill in the first pass; bounded tools; ASCII markers; depth 3; own skill; role block only | 2026-09-19 |
| K8 | PROPOSED, UNBLOCKED (2.1.278; F10 measured −71K/dispatch) | F10 | 2026-09-20 |

Rulings already given during the evaluation session (not on proposals, on framing): the substrate must satisfy Fable AND Opus (2026-09-05, → D1–D3); roles: Fable/Opus develop code with steering, Haiku categorises in fan-outs, Sonnet unused (2026-09-05, → D4); the goal is improvement by root-cause analysis and procedural checklists, not only size (2026-09-05, → Group B).

---

# Part IV — Open questions for discussion

1. [RULED 2026-09-19: docs] 1. **Evidence home** — `.claude/evidence/` (recommended: greppable, sub-agent-readable, outside the recursive `rules/` glob) vs. keeping cases in the rule files inside HTML comments (stripped from CLAUDE.md at injection per the doc; unverified for `rules/`; a Read then shows both).
2. [RULED 2026-09-19: always-on for now] 2. **plan-authoring core: always-on** (recommended — §2/§4 govern issues, commit messages and briefs) **vs. path-scoped** to `.claude/plans/**` (saves ~7K per non-plan dispatch).
3. [RULED 2026-09-19: the rule stays, the lesson goes] 3. **Which twin survives** where a lesson and a rule clause state one mechanism (18 candidates in `lessons.index.md` appendix i). Recommendation: the RULE clause is the always-on spelling; the lesson keeps only its `twin:` pointer.
4. [RULED 2026-09-19: recall test, K5] 4. **Run the optional A/B** (≈2 × 160K tokens per model) or rely on G5 alone?
5. [RULED 2026-09-19: keep IDs] 5. **Reorganise by root while keeping § numbers as IDs**, or renumber and rewrite every citation (82 rows, agent memories, issues)? Recommendation: keep IDs, add the root layer above them.
6. [RULED 2026-09-19: every key agent; VII.2] 6. **Which agents may spawn**, and the depth cap: recommendation `qa` + `numerics-investigator` only, depth 2. Alternative: depth 1 (nesting off) and everything through the main agent as today.
7. [RULED 2026-09-19: do not merge] 7. **Merge `coding-elegance` and `nexus-elegance`** into one root → tool → judgment skill?
8. [RULED 2026-09-19: rule + skill, K1, first pass] 8. **The shared instrument doctrine (B4): a rule** (always-on, ~2K, every model every dispatch) **or a skill** (loaded by the protocol and preloaded by the V&V agents)?
9. [RULED 2026-09-19: hook] 9. **claim_lint as a hook** (PreToolUse on `gh issue create` / `git commit`, warning not deny) or as a tool the checklist names?
10. [RULED 2026-09-19: ASCII markers, K4] 10. **Glyph policy**: retire ⭐/⚠ entirely, or keep a two-symbol vocabulary (⛔ live refutation, ✅ landed)?
11. [RULED 2026-09-19: structure first] 11. **Sequencing**: structure first (P1–P4, quick, bit-identical in content) then roots (P7), or roots first so the cores are written once in their final shape? Recommendation: structure first — the cores are drafted, the acceptance instruments are ready, and the root reorganisation is a content rewrite that benefits from the smaller files being in place.
12. [SETTLED 2026-09-19: inherited, 0 tool calls — VII.3] 12. **MEMORY.md in sub-agents**: it is inherited (F4) — slim it harder than the draft (8.2 KB), or accept 2.3K per dispatch?
13. [RULED 2026-09-19: K6] 13. **The surprise-log loop under the new shape**: rows append to `evidence/plan-authoring-surprise-log.md` tagged by root (B6) — who decides when a recurring root graduates to a tool, and what is the threshold (recommendation: a root with ≥3 rows and a mechanical check gets a tool)?
14. [RULED 2026-09-19: VII.2] 14. **Role assignment (H1):** Key = test-architect, method-implementer, numerics-investigator, archivist, qa; Support = explorer, literature-researcher, cross-domain-attacker, haiku categorisers; elegance-enforcer review-only without `Agent`. Agree, or should elegance-enforcer be able to call explorer for blast radius?
15. [RULED 2026-09-19: no restriction; parent review is independent] 15. **The independence invariant (H2.1)** — "authors never dispatch their own reviewers" as a hard rule, enforced by not granting `Agent`-to-reviewer in the implementer's brief. Is there a case where an implementer legitimately calls qa mid-build (e.g. a mutation battery on a gate it just wrote)?
16. [RULED 2026-09-19: a rule for now] 16. **Encoding of the workflows (H4):** an always-on rule (~2K, every dispatch) vs. a skill loaded at session start and preloaded by Key agents only. The Orchestrator needs it every session; a Support agent never does.
17. ✅ RULED 2026-09-08 (docs; see I6 for the MD/RST split). ~~**Evidence home, revisited:**~~ `docs/development/` (RST, Sphinx-gated, Nexus-indexed, generated indexes — recommended) vs. `.claude/evidence/` (Markdown, no gates). If docs: do skill `references/` still exist for skill-MECHANICS text, with docs holding only evidence and worked cases?
18. [RULED 2026-09-19: directives; VII.4] 18. **Structured directives** (`.. lesson::`, `.. anti-pattern::`, `.. surprise::` with `:check:` / `:tell:` / `:root:` fields, Nexus-indexable like `.. error-entry::`) vs. plain labelled sections. Directives give the generator structured fields and give Nexus nodes to query; cost: extension code (in `sphinxcontrib-nexus` or a local `docs/_ext`).
19. ✅ RULED 2026-09-08 (plans stay; archive is the staging spot before exclusion; close-out records move to the docs per I6). ~~**Plans stay Markdown in `.claude/plans/`**~~ (transient), with the compaction record feeding `docs/development/` at close-out — or do active plans move to the docs too?
20. [RULED 2026-09-19 — reversed: docs → `.claude` only; see K7] ~~**Include or copy the harness-loaded cores into the docs?**~~ `include :parser:` keeps one source but makes the docs build depend on `.claude/` (a `-W` failure in a rule core breaks the docs); a generated copy is a twin. Recommendation: include, and let the docs build be the rule cores' syntax gate.
21. [RULED 2026-09-19 — data in `.nexus/ontology.toml`, the one generic directive in nexus; VII.4] ~~**Where do the typed directives live**~~ — a local `docs/_ext/devblocks.py` (fast to iterate, project-only) or `sphinxcontrib-nexus` (already owns `error-entry` and the graph-node pattern; reusable across projects; a nexus release per change)? Recommendation: prototype locally, promote to nexus once the fields settle.

---

# Part VI — Related GitHub issues, and what they change in the proposals (swept 2026-09-08)

Sweep: `gh issue list` open (200) + closed (300) on ORPHEUS, open (all) + closed (200) on `deOliveira-R/sphinxcontrib-nexus`, keyword-filtered; the six directly relevant ORPHEUS issues read in full.

## VI.1 ORPHEUS — issues that PRE-DATE and CONSTRAIN this plan

| issue | state | what it rules or measures | effect on this plan |
|---|---|---|---|
| **#308** Skill layer v2: dynamic `!` injection from the corpus + record-generated docs | OPEN (2026-07-23, **user ruling**); was blocked on #231 task #10 — **#231 CLOSED 2026-07-26**, so it is unblocked | *"a skill can carry queries against the corpus instead of copies of it — residency without drift."* Four injection channels: (1) `.inc.rst` fragments included by the page AND `!cat`-ed by the skill; (2) build-time section export by anchor into `docs/_generated/skill_fragments/`, one `_GENERATORS` row; (3) point-of-use export markers `.. skill-export:: <skill>` (a container directive); (4) live-state injection. Pilot: vv-principles → a manifest | **I2 and I7 are #308, re-derived.** I2's generator = channel 2; I7's directives = channel 3 with typed fields; the eager `error_index.md` is channel 4 already shipped. ⟹ Group I is re-labelled *"execute #308 with the root-cause content design"*; do not mint a parallel mechanism. The `_GENERATORS` table + `_make_regenerator` hook fold (`817a2b8d`) exist in `docs/conf.py` today (9 generated files) |
| **#329** Rules floor: give retirement doctrine a corpus home, then thin `coding-standards.md` to a pointer (do NOT split into a skill) | OPEN (2026-08-03) | measured the floor at 30.8 KB; retirement = 61 % of coding-standards; names `docs/architecture/` (sibling of `layering.rst`) or `docs/development.rst` as the home; *"do this when #308 runs, as a second instance of the same pattern"* | **A1/I1 for coding-standards is #329.** Its "do NOT split into a skill" ruling also applies to A3: a skill `references/` folder is acceptable only for skill-MECHANICS text; evidence goes to the corpus (consistent with the 2026-09-08 ruling) |
| **#307** parity guard for the failure-mode registry dual-homed in `principles.rst` ⇄ vv-principles | OPEN (2026-07-23) | the skill and the page are governed duplication with a prose-only sync rule | **Made moot by I2**: a generated skill core cannot drift from the page. Interim: the count-parity gate it asks for, until generation lands |
| **#332** `principles.rst` redirect list is behind the skill | OPEN (2026-08-06: 9 behind, #13–#21) | `[M]` 2026-09-08: the page still has **12** bullets; the skill is at **#36 — now 24 behind** | the drift has doubled since filing; the measurement belongs in #332 as a comment. I2 closes it structurally; the cheap gate it proposes (highest skill number vs bullet count) is the interim |
| **#380** a raw PATH in a doc is a present-tense claim no gate checks (30 stale in the error catalogue) | OPEN (2026-08-18) | paths are invisible to `-W`, to `check_docstring_xrefs.py` and to `dead_references` | applies directly to **I6**: the theory pages cite `.claude/plans/archive/*.md` by bare path today; moving those files into `docs/development/archive/` under a toctree turns the path into a `:doc:` the build validates. The path-existence gate #380 proposes is a C-group tool |
| **#422** / **#379** nested inline markup renders literal backticks (46 + 32 sites; silent at every severity) | OPEN | RST cannot nest `**bold with ``literal`` inside**` | **changes the I6 criterion**: the evidence prose (surprise log, founding cases, lessons) is exactly this markup, densely. Writing it "fresh in RST" would reproduce #422 at scale. ⟹ third leg of the criterion: prose with nested inline markup is MD (CommonMark nests emphasis and code spans), so the evidence pages are MyST Markdown with typed blocks as `:::{evidence}` fences — I6's table is amended below |
| **#455** a `catches(ERR-NNN)` marker the defect does not redden is a phantom catcher | OPEN | | an instance of X1 (the instrument must be able to fail) at the test-harness tier; a candidate row for the shared instrument doctrine (B4) |
| **#334** file-level `verifies` lists over-credit equation coverage | OPEN | | X2 at the harness tier (population of a claim) |
| #148 method-implementer: skill follow-ups | OPEN | | touched by A3 (numerical-bug-signatures split) and H5 |
| #302 silent dead py-domain xrefs (no autodoc target) | OPEN | | the docs-side twin of the retirement audit's third search; unchanged by this plan |
| #147 add the method-implementer agent | CLOSED 2026-05-02 | | precedent for agent design work; W1's build phase |

## VI.2 sphinxcontrib-nexus — relevant open items

| issue | relevance |
|---|---|
| #90 (CLOSED 2026-08-18) directive-misuse warnings logged through stdlib logging, so `-W` did not gate them | **I7's validation promise depends on this being fixed** — custom directives must raise Sphinx warnings; verify the fix covers a locally registered directive before promising "`-W` validates every typed block" |
| #87 authored edges need provenance and a dangling report | the `evidence → instances → rule` edges I7 mints are authored edges; they inherit this gap |
| #67 neighbour payloads ~5× larger than their information | the same shape as the briefing's size (F1: ≈5.5K with a 1.5K id-grammar block); a nexus-side item |
| #65 nexus-exploring/reference.md says 29 tools, there are 40 | skill drift in the nexus-shipped skills — the same defect class as #332, on the other repo |
| #39 dead-reference resolver false-DEAD classes | the instrument I3 relies on; know its two blind spots |
| #20 flow-graph directive | precedent for directive + graph work in nexus (Q21) |
| **not tracked anywhere**: the eager loading of 44 tool schemas while `session_briefing`'s `preload_hint` states they are deferred (F1, A8a) | file on nexus |

## VI.3 Amendment to I6 from #422

| durable content | format |
|---|---|
| founding cases, surprise log, lesson bodies — dense nested inline markup | **MD (MyST)**, typed blocks as colon fences; conversion from the current Markdown is a MOVE, not a rewrite (no #422 exposure, no backtick-conservation risk) |
| the vv / elegance catalogue entries (mechanism + case) — currently Markdown, same markup density | **MD (MyST)** for the same reason; the page-level doctrine (`principles.rst`) stays RST and is generated-toward, not hand-synced |

## Group J — Issue hygiene this plan owes (Cardinal Rule 4)

**J1. Comment on #308**: this plan's Group I is its execution design; link the decision document; note #231 is closed so the blocker is gone. **J2. Comment on #332** with the 2026-09-08 measurement (24 behind). **J3. Comment on #329** that the corpus home is `docs/development/` per the 2026-09-08 ruling. **J4. File on nexus**: eager tool schemas vs `preload_hint`. **J5. File on ORPHEUS** (after the rulings): one umbrella issue for the substrate restructure carrying this document's proposal IDs, labelled `module:docs`, so a fresh session can pick it up from the issue. — Status: PROPOSED (J1–J3 are comments on existing issues and could be posted now; J4–J5 wait for the rulings).

---

# Part VII — Rulings round 1 (2026-09-19) and what they produced

## VII.1 Rulings banked (verbatim where the user's words carry the criterion)

| Q | ruling | consequence |
|---|---|---|
| 1 | evidence is durable ⟹ docs; RST vs MD per the I6 criterion (mine to apply) | I1/I6 stand |
| 2 | "for now, we keep plan_authoring in context"; revisit if a better pattern or harness option appears | A1 keeps plan-authoring always-on |
| 3 | "A rule is a lesson that has been distilled and lifted. The rule stays, the lesson goes. A twin is wasted context." | A4: the 18 absorbed lessons are RETIRED from the index (bodies to the docs evidence page); no `twin:` pointers survive |
| 5 | keep § numbers as stable IDs, add the root layer above | B1 as proposed |
| 6, 14 | "All key sub-agents should be capable of spawning support agents" — a method-implementer spawning a numerics-investigator is legitimate (it keeps the parent's context clean); **elegance-enforcer is a KEY agent** (final parallel review with qa); both may spawn support | H1/H5 amended (VII.2) |
| 7 | do NOT merge `coding-elegance` (what elegance is) with `nexus-elegance` (how to use the tool to investigate it) | B3 amended: root → tool → judgment stays a two-file pair; Q7 closed |
| 9 | claim lint as a hook | C1 → PreToolUse hook on `gh issue create/comment` and `git commit`, warn not deny |
| 10 | glyphs were never the user's; acceptable iff the meaning is consistent and unambiguous | K4 (marker vocabulary) |
| 11 | sequencing accepted; "we will do multiple passes until it looks good" | G8 as proposed; expect ≥2 passes |
| 12 | first eliminate "did it search for it?" | settled, VII.3 |
| 13 | check the tool count does not grow without bound; teach HOW TO DO IT RIGHT (narrow); teach how-not-to only for a recurring unpleasant pattern | K6 |
| 15 | no restriction on a sub-agent spawning qa: "even if the qa the method-implementer spawned says all good, the parent spawns a different independent qa anyway" — a spawn restriction would only matter if a sub-agent could rule on its own work | H2.1 amended: the invariant is *the parent's review is independent*, not *authors may not spawn reviewers* |
| 16 | workflows as a RULE for now; review after the first reorganisation | H4(a) as proposed |
| 18 | structured directives; extension code local; check the nexus ontology mechanism | VII.4 |
| 20 | **"The flow should be from docs to .claude, never the opposite … we are preparing a harness-independent resource, such that your brain is not tied to the claude code harness"** | REVERSES my Q20 recommendation: nothing in `docs/` includes `.claude/`; `.claude/` is GENERATED from the docs — K7 |
| 21 | check Nexus's node/edge extension; we own sphinxcontrib-nexus | VII.4 |
| extra | an ARTICULATION rule (standards of articulation + no mannered prose); hottest in plan writing and documentation writing | K2 |
| extra 2 | CLAUDE.md: the Cardinal Rules become rules; CLAUDE.md becomes an ON-BOARDING file (architecture + direction), crafted together | K3 |

## VII.2 Amendments to Group H

- **H1 roles (amended).** Key agents: test-architect · method-implementer · numerics-investigator · archivist · **qa · elegance-enforcer**. Support: explorer · literature-researcher · cross-domain-attacker · haiku categorisers. Every KEY agent holds `Agent`; every SUPPORT agent does not. Review is a PHASE (qa ∥ elegance-enforcer, dispatched by the parent on the artefact), not a role that forbids spawning.
- **H2.1 (amended).** *The parent's review is independent of the child's.* A child may spawn any agent it needs, including qa; the parent still dispatches its own qa and elegance-enforcer on the result. Nothing is forbidden; independence is guaranteed at the parent, where it is cheap to guarantee.
- **H2.2 (amended 2026-09-19, Q22) depth ≤ 3 is the CEILING, 2 the norm** — [RULED] depth 3 so an implementer's numerics-investigator keeps its supports; the harness default is 3, so nothing is written. Superseded text follows for the record: ~~depth ≤ 2~~ — orchestrator → key → support — enforced by `CLAUDE_CODE_MAX_SUBAGENT_SPAWN_DEPTH=2` and support allowlists; a key agent's spawn of a KEY agent (implementer → numerics-investigator) is then the deepest legal chain, and the investigator's own supports are withheld by the cap. ⚠ Open: does depth 2 leave the investigator able to spawn explorer? By the doc, a depth-2 agent has no `Agent`. Options: depth 3 with the cost accepted for that one chain, or depth 2 and the investigator does its own exploring. → Q22.
- **H5 (amended).** Add `Agent` to test-architect, method-implementer, archivist, elegance-enforcer; keep on numerics-investigator, qa; keep OFF explorer, literature-researcher, cross-domain-attacker. The Opus delegation sentence + a return-payload cap on every key agent.

## VII.3 Q12 settled — MEMORY.md IS inherited by sub-agents in this harness

`[M]` both probes' transcripts contain **0 `tool_use` blocks** (grep over the JSONL, not read); the second probe quoted, verbatim, the 80 characters FOLLOWING the phrase I asked about — text that was not in my prompt — and the first Active-work bullet. No search happened; the text was in its context. ⟹ MEMORY.md is in the always-on budget for every dispatch (2.3K after A5). Decision needed only on how hard to slim it: the A5 draft (8.2 KB) is under half the cap; the four disciplines plus "hooks ≤ 15 words" are the mechanism, not a one-off cut.

## VII.4 Q18/Q21 — what Nexus can extend as DATA, and the one piece of CODE it needs

`[M]` read `sphinxcontrib/nexus/ontology.py`, `ontology.toml`, `directives.py`, `.nexus/ontology.toml`:
- **Node and edge TYPES are data.** A project declares `[node.<name>]` / `[edge.<name>]` in `.nexus/ontology.toml` (monotone widening; base entries cannot be redefined); `Ontology.node_types` includes project-declared types "which is the entire point of the extension tier"; the graph stores `type` as a string. ORPHEUS already uses the mechanism (`[extend.edge.implements] domain = ["data"]`, 2026-08-18).
- **Minting is code.** The docstring is explicit: *"Declaring `exercises: test → code` costs three lines; populating it still needs an extractor."* `apply_declared_nodes` handles exactly one directive kind (`error-entry`, hard-coded); a directive enqueues `{"kind","id","title","docname","lineno"}` on `env.nexus_pending_edges` and the mint runs at build-finished.
- ⟹ **The right split:** ORPHEUS's `.nexus/ontology.toml` declares `dev_rule`, `dev_evidence`, `dev_lesson`, `dev_antipattern`, `dev_rationale` (+ edges `instances: dev_evidence → dev_rule`, `rationale_for: dev_rationale → equation`) — data, no release. **sphinxcontrib-nexus gains ONE generic declaring directive** (`.. nexus-node:: <type> <id>` with `option_spec` derived from the type's declared `attributes`, rendering as an admonition; or a factory registering one directive per ontology-declared type) and a generalised `apply_declared_nodes` keyed on the ontology instead of on `"error-entry"`. Local ORPHEUS extension code then shrinks to rendering preferences at most. This is a nexus issue to file (J6), and it also retires the hard-coded `error-entry` special case (elegance C2: discriminate once, on the ontology).
- Q21 answer: the directive machinery lives in nexus (we own it; it already owns the pattern); the vocabulary lives in ORPHEUS's ontology file.

## VII.5 Rulings round 2 (2026-09-19, later)

| item | ruling | consequence |
|---|---|---|
| Q4 | recall test accepted | G6 replaced by K5 |
| Q8 | rule + skill accepted, "in our first implementation" | K1 is in scope for the first pass |
| Q13 | agreed | K6 governs Group C |
| Q10 | ASCII markers accepted | K4; glyphs retired from all generated surfaces |
| Q22 | depth 3 | `CLAUDE_CODE_MAX_SUBAGENT_SPAWN_DEPTH=3` (the default — so no setting is written; H2.2 amended to say depth 3 is the ceiling, 2 the norm) |
| Q23 | its own skill | `instrument-doctrine` skill, preloaded by qa, test-architect, numerics-investigator, archivist; vv-principles cites it |
| Q24 | role block only | K7: AGENT.md headers are harness-specific (e.g. `omitClaudeMd` exists only from 2.1.271); only the role block is generated from `docs/development/workflows` |

**All 24 open questions are now ruled or settled.** The plan is executable from P1. The user's remark that sub-agents can now omit CLAUDE.md produced F10 and K8.

## Group K — proposals from this round

**K1. Instrument doctrine = a RULE (what) + a SKILL (how).** (Q8)
- *Rule* `.claude/rules/instrument-doctrine.md`, always-on, ≤ 800 tokens, four statements each with its one check: **X1** an instrument is evidence only if some realizable state changes its reading — before citing a gate, metric, canary or fix, name the state that reddens it and show it ran; **X2** every claim carries its population and its instrument — write `k of N <predicate>`, the command, the fixture, the draw/protocol; **X3** prose is not enforcement — a docstring, marker, label, table header or plan sentence asserts nothing the code or gate does not; assert the structure the prose names; **X4** one definition per quantity — two spellings of one thing agree tautologically; find the shared upstream before calling agreement evidence.
- *Skill* (`instrument-doctrine`, or the first section of `vv-principles` core): the procedures — positive control (in-class mutation, per-arm table); stabiliser enumeration for a functional; activation counting for a canary; denominator + predicate + validated filter for a census, the two-filter completeness check; draw-stable statistics; the α-normalised AST check for "independent" implementations. Loaded by the protocol batch and preloaded by qa, test-architect, numerics-investigator, archivist.
- The cores of plan-authoring, vv-principles, coding-elegance, coding-standards cite X1–X4 by ID and keep only their domain-specific instances — the cross-file duplication in F5 collapses to one page. Why rule+skill: the four statements apply to every artefact an agent writes (issues, commits, briefs), so they must be always-on; the procedures are needed only when building or judging evidence.

**K2. Articulation rule** `.claude/rules/articulation.md` (promotes the user's `feedback_articulation_lossless_disassembly` memory). Content: (i) *articulate = disassemble a concept so a fresh reader reassembles it losslessly; the LOSS is the measure* — for a plan the loss is a surprise, for a doc a dead reference or a reader's wrong reconstruction; (ii) *no mannered prose* — the Fable 5.1 page's definition ("say what you mean; when a literal phrase is available, use it"); metaphor drags in connotations you did not choose; it also costs context; (iii) straight to the point is not clipped — one idea per sentence, the reader's vocabulary not the session's, names spelled out once; (iv) the two hot cases: a plan is a message to yourself or a sub-agent after context is gone (plan-authoring governs its claims; this rule governs its prose), documentation is the brain (Cardinal Rule 3); (v) the check — the vendor's golden rule: *show it to a colleague with minimal context; if they would be confused, so will the model*; plus the plan-side check: every named object is either defined in the text or linked. ≤ 600 tokens. Status: PROPOSED.

**K3. CLAUDE.md → rules + an on-boarding file.** (extra 2) Proposed split: Cardinal Rules 1–5 → `rules/correctness.md` (1), `rules/architecture.md` (2, with the coding-elegance pointer), `rules/sphinx-brain.md` (3, with the archivist/explorer triggers), `rules/issues-are-the-log.md` (4), `rules/delegation.md` (5 merges into the existing file + the workflows rule H4); "Working Principles" → articulation (K2) + coding-standards; "Git Workflow" → docs/development/git_workflow + a PreToolUse hook for the two mechanical items; "Environment resolution" → `rules/environment.md`. CLAUDE.md becomes the on-boarding page, crafted with the user: what ORPHEUS is and for whom; the architecture map (the six core folders, the layer contract, the operator-algebra spine, where the docs live); the direction of development (SN first as the vanguard, the sharpening order, what is deliberately not being invested in); how a session starts (the protocol pointer); where the rules and workflows are (a generated index). Target ≤ 1.5K tokens. Status: PROPOSED — skeleton only, content is a joint session.

**K4. Marker vocabulary — one meaning, one spelling, typable.** (Q10) The user types on a terminal and never used glyphs; I introduced ⭐/⚠/⛔/⟹. Proposal: canonical markers are ASCII brackets already in use — `[M]` measured (+ command/date), `[R]` reasoned (not yet measured), `[REFUTED <date>]`, `[LANDED <hash>]`, `[REMEDIED <date> @<hash>]`, `[HYPOTHESIS]` for a means proposed before investigation. Glyphs are RETIRED from all cores and generated surfaces (⭐ "important" carries no information; ⚠ and ⛔ overlap; ⟹ is an arrow chain the Fable 5.1 page asks to avoid). Rendering in the docs may map markers to admonition classes. A one-line table in the articulation rule fixes the meanings. Status: PROPOSED.

**K5. Replace the n=1 A/B with a RECALL test against a known finding set.** (Q4 explained) The A/B proposed in G6 runs the same review twice — full skill vs core — and compares finding sets; with one run per arm the difference is indistinguishable from run-to-run variance, and it costs ≈2 × 160K per model. The alternative: pick a PAST review whose findings are recorded and were confirmed by the fix that followed (e.g. the #426 review round, `f52877db`, or any qa report in `scratch/` with a landed fix), check out the pre-fix commit, dispatch qa with the CORE skill on that diff, and count how many of the known findings it reproduces — `k of N`, one run per model, against ground truth instead of against another draw. Recommendation: K5, on both models (Opus as the qa preload, Fable as a main-session review), and keep G5 (surprises per campaign, by root) as the ongoing detector. Reason: a recall test has a denominator; an A/B at n=1 does not.

**K6. Tools encode the POSITIVE procedure; anti-patterns attach as tells. Bounded growth.** (Q13) A tool per recurring FAILURE would grow with the failure catalogue (unbounded: how-not-to is the wide space). Instead every tool implements one item of a floor checklist — the narrow "how to do it right" — so the count is bounded by the checklist: plan floor 11 items, verification floor 7, elegance 14; of these `[R]` ≈ 8 are mechanical (claim carries population → lint; blast set by contract → census; date collision → check; canary activation → counter; anchor/reference validity → the docs build; path existence → #380's gate; marker set → the harness audit; retirement residue → `dead_references`). A recurring negative pattern (≥ 3 surprise-log rows on one root) is recorded as a TELL on the existing checklist item and, if mechanical, becomes a CASE in that item's tool — never a new tool. Expected count: ≈ 8 tools, then flat. Status: PROPOSED.

**K7. `.claude/` is generated from `docs/development/`; the docs never depend on `.claude/`.** (Q20 reversed) The source of every rule core, skill body and agent role block lives in `docs/development/` (MyST Markdown where the content is prose with typed fences, RST where it owns equations/directives — the I6 criterion applied to SOURCES); a generator (`tools/docs/generate_harness.py`, one `_GENERATORS` row at `build-finished`, `--check` in CI) emits `.claude/rules/*.md`, each `SKILL.md` body under its frontmatter, and the per-agent role block, stamped `GENERATED — edit docs/development/…`. The generated files are COMMITTED (the harness reads them without a build), exactly like `error_index.md` today. Where the source is already Markdown the "generation" is a stamped copy plus a token-budget check (a core over its budget fails `--check`), so the source IS the artefact and typed fences (`:::{evidence}`) survive verbatim as the content-type separators the vendor doc asks for. Evidence pages have NO `.claude/` copy — agents read `docs/development/*.md` on demand. Any other harness (an API caller) reads the same `docs/development/` tree. I6's "include `.claude` into the docs" row is REFUTED by this ruling and struck. Status: PROPOSED (mechanism); the premise is RULED.

**K8. `omitClaudeMd: true` on every SUPPORT agent — UNBLOCKED 2026-09-20 (F10 measured: −71K per dispatch).** Support agents (explorer, literature-researcher, cross-domain-attacker, haiku categorisers) take everything from the brief by definition (H1); they do not need the rules, and the rules are the bulk of their cost. Key agents keep the block. Prerequisite met: Claude Code 2.1.278 installed; F10 shows the switch drops CLAUDE.md + all rules + MEMORY.md (≈71K). The brief template for Support agents (H2.4) must then carry the two or three rules a Support agent still needs (e.g. the ugrep hazard for explorer; the local-folder-first rule for literature-researcher) — a Support brief is the only rule surface those agents see. Status: PROPOSED — ready; recommended for P1 alongside the rules split, since it is a four-line frontmatter change on four agents plus the Support brief template (H2.4).

**J6. File on sphinxcontrib-nexus:** a generic ontology-driven declaring directive + `apply_declared_nodes` keyed on the ontology (VII.4). Status: PROPOSED.

22. [RULED 2026-09-19: depth 3] 22. **Spawn depth 2 or 3?** With depth 2, an implementer that spawns a numerics-investigator leaves the investigator without `Agent` (no explorer, no literature pull). Options: depth 3 and accept the cost for that chain; or depth 2 and the investigator explores itself. Recommendation: depth 3, since the cap is a ceiling not a norm, the investigator's supports are precisely where a fresh context pays, and the always-on block shrinks 64K → 24K under this plan.
23. [RULED 2026-09-19: its own skill] 23. **Where does the instrument doctrine's SKILL half live** — its own skill (`instrument-doctrine`, preloaded by the four V&V agents) or the first section of the vv-principles core? Recommendation: own skill; vv-principles then cites it and keeps only the V-specific instances, which is the same single-source move as K1's rule half.
24. [RULED 2026-09-19: role block only; headers are harness-specific] 24. **K7 details**: are AGENT.md files generated in full, or only their role block? Recommendation: role block only in the first pass (the rest of an AGENT.md is harness-specific: tools, model, memory).

---

# Part VIII — P1 execution record (2026-09-20, branch `docs/development-substrate`)

**What landed** (one commit; every path listed in the commit body):
- `docs/development/` — the harness-independent source: `index.rst`, `git_workflow.rst` (moved), `harness.md`, `workflows.md`, `lessons.md` (44-line index), `rules/` (plan-authoring, coding-standards, instrument-doctrine, articulation, workflows), `skills/` (vv-principles, coding-elegance, instrument-doctrine), `agents/` (9 role blocks + index), `evidence/` (plan-authoring: 94 surprise-log rows + 25 founding cases; coding-standards: 18 cases; lessons: 64 bodies; vv-anti-patterns: 36; test-design-modes: 6 + 2 worked cases; coding-elegance: 51 sections). All evidence moved VERBATIM by script (each agent asserted line-set identity).
- `tools/docs/generate_harness.py` + `harness_manifest.toml` (18 targets; `--check`: drift, budget, dead heading link; a `_GENERATORS` row in `docs/conf.py`); `myst-parser>=5` in the docs extra; `conf.py`: `myst_parser`, colon fences, `myst_heading_anchors = 4`.
- `.claude/` GENERATED: 5 rule cores, 3 skill cores, `lessons.md`, 9 role blocks. Retired: `subagent-handoff-protocol` skill (its return contract is `docs/development/workflows.md`); the 153 KB `lessons.md` (bodies in the docs); glyphs on every core (K4).
- Agents: `Agent` on the six Key agents; `omitClaudeMd: true` on the three Support agents; `instrument-doctrine` preloaded by qa, test-architect, numerics-investigator, archivist; stale handoff pointers repointed (CLAUDE.md, delegation.md, method-implementer). `session-start.txt` names the new batch. `settings.json`: the inert `Write(.claude/agents/**)` rule removed. MEMORY.md re-indexed 17.6 → 8.8 KB (backup `scratch/_harness_eval/backup/MEMORY.md.20260920`).
- Corpus repairs the move forced: 11 `:doc:` citations of the old `development` page repointed to what they meant (vv-principles / coding-elegance / the error catalogue / a literal skill name); 2 relative `:doc:` links in the moved page made absolute; one H2 inserted in an evidence page (MyST header-level rule); one 12 373-char paragraph soft-wrapped (docutils line limit).

**Measured after generation** (`[M]` chars/3.6 over the generated files):

| surface | before | after |
|---|---:|---:|
| plan-authoring rule | 40.3K | 8.9K (110 mechanisms) |
| coding-standards rule | 9.6K | 4.9K |
| new rules: instrument-doctrine / articulation / workflows | — | 1.0K / 0.7K / 1.2K |
| untouched rules: process-discipline / nexus-tools / delegation / CLAUDE.md | 3.3K / 2.5K / 1.3K / 2.7K | unchanged (second pass: K3) |
| **always-on block** (rules + CLAUDE.md, ex. path-scoped vv-testing) | **≈64.6K** | **≈26.8K** (+ MEMORY.md 2.5K) |
| vv-principles skill | 42.5K (+2.6K injected) | 9.4K (36 + 6 items; injection kept) |
| coding-elegance skill | 21.6K | 8.7K |
| instrument-doctrine skill | — | 1.9K |
| lessons | 40.9K | 3.1K |
| **session-start batch** (+ briefing ≈5.5K) | **≈125K** | **≈28.7K** |
| Support-agent dispatch floor (haiku, zero tools) | 107.5K | 36.5K (`omitClaudeMd`) |

**Gates:** `generate_harness --check` 0 problems / 0 drift; harness-facing tests (`test_error_catalogue_reconciles`, `test_docstring_xrefs`, `test_elegance_debt_is_tagged`, `test_layer_imports`) 421 passed under `-O`; `sphinx -W`: first run 31 warnings in the three classes above, all fixed; second run: **0 warnings, rc=0**; Nexus `dead_references` on the rebuilt graph: 0 dead of 66 checked.

**Left for the second pass (P1b), by ruling "multiple passes":** K3 (CLAUDE.md → rules + on-boarding page, a joint session); the untouched rules under generation; `principles.rst` redirect list → the evidence page (#332); root-layer reorganisation (B1/B2, § IDs kept); typed directives (I7, after the nexus J6 issue); the Support brief template into the agents that dispatch them; `[[lessons-Lnn]]` wiki-links in the elegance core wired to the lessons anchors.

---

# Part V — Evidence appendix

**Commands and probes (reproducible).** Sizes: `wc -c` over the surfaces in F1; section splits by `re.split(r'(?m)^(?=#{1,3} )')`; surprise-log rows by `^\| 2026`; repeat/no-new-clause by regex on the clause column; growth by `git log --format=%ad --date=format:%Y-%m -- <file>` and `git show <hash>:<file> | wc -c` at month-end commits; emphatic tokens `grep -oE '\b(CRITICAL|NEVER|MUST|ALWAYS|DO NOT|IMPORTANT)\b'`; glyphs `grep -oE '⭐|⛔|⚠'`; agent preloads from `skills:` frontmatter; models from `model:` frontmatter and `git log -400 --format=%b | grep -oE 'Co-Authored-By: Claude [A-Za-z]+ [0-9.]+'`. Probes: (i) haiku `general-purpose`, ten YES/NO strings, zero tools — 113 006 tokens, all five rule files + MEMORY.md present, 43 Nexus tools visible; (ii) haiku, verbatim quote of MEMORY.md's first Active-work bullet — reproduced; (iii) sonnet → haiku → haiku nested dispatch — `PONG from depth 3`, depth-3 tool list without `Agent`.

**Pages read 2026-09-08:** myst-parser.readthedocs.io (roles-and-directives: fences, roles, `(label)=`, `eval-rst`); docutils.sourceforge.io/docs/ref/rst/directives.html (`include` `:parser:` — docutils ≥ 0.17, provisional; `admonition` options `class`/`name`; `container`, `class`, `topic`, `sidebar`); sphinx-doc.org directives (`only` = build tags via `-t`, "control only content of document").

**Vendor pages read (2026-09-05).** code.claude.com/docs/en/memory (200-line guidance; recursive `rules/`; path-scoped on Read; HTML comments stripped; `InstructionsLoaded`; MEMORY.md cap; "auto memory not loaded into subagents" — contradicted by probes); code.claude.com/docs/en/sub-agents (3-layer nesting default; depth/concurrency env; `tools:` allowlist withholds `Agent`; `SendMessage`; resume by name); platform.claude.com …/claude-prompting-best-practices (XML for content types; add the WHY; `<example>` tags; long-context ordering); …/prompting-claude-fable-5 (brief over enumeration; too-prescriptive skills degrade output; one lesson per file; explicit verification cadence; no reasoning-echo instructions); …/prompting-claude-fable-5-1 (writing density; formatting; finish the whole task; compaction summary contents; scope of changes; lead agent keeps working); …/prompting-claude-opus-5 (longer responses and files; remove verification/re-check instructions; literal pre-filters; delegation caps; correction narration); …/agent-skills/best-practices (concise; degrees of freedom; SKILL.md ≤ 500 lines; references one level deep; TOC > 100 lines; no time-sensitive info; evaluations first; test with all models).

**Draft inventory** (`scratch/_harness_eval/`): `plan-authoring.core.md`, `coding-standards.core.md`, `MEMORY.index.core.md`, `vv-principles.core.md`, `coding-elegance.core.md`, `subagent-handoff.core.md`, `lessons.index.md`, `session-start.txt.proposed`, `plan-checklist.md`, `worklog_2026-09-05.md` (the chronological record with the full root-assignment tables §9–§10 and the Fable/Opus reconciliation §8b).

**Growth and line counts.** plan-authoring 1 022 lines (143 chars/line); lessons.md 2 475; vv-principles 1 788; subagent-handoff 733; coding-elegance 681; coding-standards 451; process-discipline 185; CLAUDE.md 182; nexus-tools 131; delegation 70; MEMORY.md 63.
