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
- *Rule* `.claude/rules/instrument-doctrine.md`, always-on, ≤ 800 tokens (`[M]` 2026-09-20 ≈1003, budget 1100: the cap was a proposal, the manifest budget is the instrument), four statements each with its one check: **X1** an instrument is evidence only if some realizable state changes its reading — before citing a gate, metric, canary or fix, name the state that reddens it and show it ran; **X2** every claim carries its population and its instrument — write `k of N <predicate>`, the command, the fixture, the draw/protocol; **X3** prose is not enforcement — a docstring, marker, label, table header or plan sentence asserts nothing the code or gate does not; assert the structure the prose names; **X4** one definition per quantity — two spellings of one thing agree tautologically; find the shared upstream before calling agreement evidence.
- *Skill* (`instrument-doctrine`, or the first section of `vv-principles` core): the procedures — positive control (in-class mutation, per-arm table); stabiliser enumeration for a functional; activation counting for a canary; denominator + predicate + validated filter for a census, the two-filter completeness check; draw-stable statistics; the α-normalised AST check for "independent" implementations. Loaded by the protocol batch and preloaded by qa, test-architect, numerics-investigator, archivist.
- The cores of plan-authoring, vv-principles, coding-elegance, coding-standards cite X1–X4 by ID and keep only their domain-specific instances — the cross-file duplication in F5 collapses to one page. Why rule+skill: the four statements apply to every artefact an agent writes (issues, commits, briefs), so they must be always-on; the procedures are needed only when building or judging evidence.

**K2. Articulation rule** `.claude/rules/articulation.md` (promotes the user's `feedback_articulation_lossless_disassembly` memory). Content: (i) *articulate = disassemble a concept so a fresh reader reassembles it losslessly; the LOSS is the measure* — for a plan the loss is a surprise, for a doc a dead reference or a reader's wrong reconstruction; (ii) *no mannered prose* — the Fable 5.1 page's definition ("say what you mean; when a literal phrase is available, use it"); metaphor drags in connotations you did not choose; it also costs context; (iii) straight to the point is not clipped — one idea per sentence, the reader's vocabulary not the session's, names spelled out once; (iv) the two hot cases: a plan is a message to yourself or a sub-agent after context is gone (plan-authoring governs its claims; this rule governs its prose), documentation is the brain (Cardinal Rule 3); (v) the check — the vendor's golden rule: *show it to a colleague with minimal context; if they would be confused, so will the model*; plus the plan-side check: every named object is either defined in the text or linked. ≤ 600 tokens (`[M]` 2026-09-20 ≈708, budget 800). Status: PROPOSED.

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

---

# ⏸ COMPACTION POINT #1 — 2026-09-20. P1 LANDED, UNMERGED; next session REVIEWS before merge

**Read this section first after compaction.** Then Part VIII (what landed, sizes, gates), then Part VII (the rulings). Do not re-litigate a ruling; every open question 1–24 is ruled.

**The order of work after compaction is the Transition and evaluation protocol (T1–T8) at the end of this section: T1 (the review charter) and T2 (the `InstructionsLoaded` hook) happen in the post-compaction session; T3 is the restart ON THE BRANCH; T4–T6 measure; T7 merges; T8 resumes the queue.**

**[LANDED 0d376e4f, second pass in the next commit] T1 and T2 are DONE — the record is the section "T1 and T2 executed" and its "second pass" subsection at the end of this compaction point. T3 (the restart on the branch) and T4 (the measurements) are DONE `[M]` 2026-09-20 — the record is the section "T4 executed": always-on block 30 655 tokens (target 25K missed by the four untouched files), Key-shaped floor 74 822 (target 75K met), Support 44 167. T5 (the planted task) is DONE `[M]` 2026-09-20 — the record is the section "T5 executed": one Support and one Key dispatch, both contracts honoured, two surprises by root (L38; §2 [M]-SCOPE), a false premise in the paragraphs caught by the qa's own probe. T6 (the recall test) is DONE `[M]` 2026-09-20 — the record is the section "T6 executed": k of 9 known findings reproduced, F-1 excluded as contaminated by the always-on rule on both substrates — old_qa 5, new_qa 6.5, old_fable 6, new_fable 5.5; one draw per cell. T7 (the merge) is DONE 2026-09-20 — the record is the section "T7 executed": #94 and #95 on sphinxcontrib-nexus, #477 the umbrella, #308 commented, ff-merged and pushed. The next step is T8: the queue in #477, on `main`.**

## Where things stand

| item | state | how to verify |
|---|---|---|
| P1 commit | `[LANDED 72006892]` on `docs/development-substrate`, 62 files, +9 881 / −7 211 | `git log -1 --stat 72006892`; `git merge-base --is-ancestor 72006892 main` is FALSE until merged |
| tree | clean except untracked `scratch/` (a `git clean` destroys the drafts, the gate scripts, the MEMORY.md backup) | `git status --porcelain \| grep -v '^?? scratch/'` → empty |
| generator | 18 targets, 0 problems, 0 drift | `.venv/bin/python -m tools.docs.generate_harness --check` |
| docs build | 0 warnings, rc 0 (second run; first had 31, all repaired) | `.venv/bin/python -m sphinx -W --keep-going -q -b html docs docs/_build/html` (≈10 min; run detached, log to `scratch/_harness_eval/`) |
| Nexus | `dead_references` 0 of 66 on the rebuilt graph | `mcp__nexus__dead_references` |
| tests | 421 passed: `tests/test_error_catalogue_reconciles.py test_docstring_xrefs.py test_elegance_debt_is_tagged.py test_layer_imports.py` under `python -O` | same command; the change touches no `orpheus/`, so the full gate is not owed by this commit |
| MEMORY.md | re-indexed 17.6 → 8.8 KB; backup `scratch/_harness_eval/backup/MEMORY.md.20260920` | outside git |
| this session's harness | still runs on the OLD rules/CLAUDE.md snapshot (session-start snapshot); the generated cores are live only from the next session; skill and AGENT.md edits ARE live now | `reference_harness_context_snapshot_timing` |

## The review charter for the post-compaction session (before merge)

The user's instruction: review the work, run tests and reviews, and merge only when satisfied with the quality; then compact again and restart. Reviews are dispatched by the parent (H2.1), in parallel where independent; every brief carries the return contract (word cap; files carry the detail; `NEEDS:`).

1. **Fidelity of the cores to the rulings** — read, do not delegate: `docs/development/rules/{instrument-doctrine,articulation,workflows}.md` (written by the main agent from K1/K2/H1–H4; nobody else has read them) and `docs/development/workflows.md`, `harness.md`, the nine `agents/*.md`. Check each against Part VII's rulings verbatim (roles, invariants, depth 3, parent-review independence, the return contract, ASCII markers). Fix in the SOURCE and regenerate.
2. **Fidelity of the distilled cores to their originals** — dispatch `qa` on `docs/development/rules/plan-authoring.md` vs `git show main:.claude/rules/plan-authoring.md`, and on `skills/vv-principles.md` vs `git show main:.claude/skills/vv-principles/SKILL.md`: does every clause / anti-pattern / mode survive with its imperative, check and tell? Are any two merged that catch different failures? The agents' own appendices (110 mechanisms; 36 + 6 items) are the claimed denominators; qa recounts. Same for `coding-standards` and `coding-elegance` at lower priority (their agents reported no post-draft deltas).
3. **Verbatim-move audit of the evidence pages** — one `general-purpose` agent with a script: for each evidence page, every non-heading, non-`**Surprise.**/**Clause.**/Clause:` line must appear in the corresponding `main` original (allowing the documented normalisations: `\|` → `|`, un-quoted `>`, the one soft-wrapped paragraph in `test-design-modes.md`, the one link-to-code-span in Mode 12). Report `k of N` lines matched per page.
4. **The generator** — dispatch `elegance-enforcer` on `tools/docs/generate_harness.py` + the manifest (a 200-line tool: link re-pointing, slugging, agent-block insertion, `--check`); and `qa` for a mutation check: break a heading in an evidence page → `--check` must report the dead link; edit a generated file by hand → `--check` must report drift; push a core over budget → must fail.
5. **The eleven repointed citations** (Part VIII) — read each in context (`git show 72006892 -- docs/theory/…`): does the new target say what the sentence claims? Two now point at the `vv-principles` skill page for "Bit-identity", "1-group degeneracy", "Structural independence" — confirm those sections exist on the page (they do at `skills/vv-principles.md` §§ "Bit-identity vs principled-equivalence", "1-group degeneracy", "The three pillars") and that `theory/verification/principles.rst` is not the better target (it is the corpus doctrine page; if so, repoint there instead).
6. **Agent frontmatter matrix** — re-print it (the script in the session worklog; or `grep -A20 '^tools:' .claude/agents/*/AGENT.md`) and confirm against H5: Agent on the six Key agents; `omitClaudeMd: true` on the three Support agents; `instrument-doctrine` preloaded by qa, test-architect, numerics-investigator, archivist; no `subagent-handoff-protocol` anywhere (`grep -rn subagent-handoff-protocol .claude docs/development` → only the History paragraph in `workflows.md`).
7. **Support-agent briefs** — the three `omitClaudeMd` agents see NO project rule. The brief template (`docs/development/workflows.md` § The brief) carries the rules; but NOTHING yet makes a Key agent include them. P1b item: a "Support briefs" paragraph in each Key agent's role block. Decide whether that must land BEFORE merge (recommendation: yes, it is ten lines per agent and the exposure is real — an explorer census brief without the ugrep hazard).
8. **Then the gates again**: generator `--check`; `sphinx -W`; `dead_references`; the four test modules; and `git diff main --stat` read in full (§ the `git add -A` rule: the staged set was explicit, 62 paths; confirm nothing under `scratch/` is in the commit: `git show --stat 72006892 \| grep scratch` → empty).
9. **Do NOT merge yet** — the merge is T7 of the Transition and evaluation protocol below, after the fresh-session measurements T4–T6 on the unmerged branch. After step 8 is green: T2 (the hook), then compact, then restart ON THE BRANCH.

## Status deltas this compaction records (Part II statuses are otherwise unchanged)

IMPLEMENTED @ `72006892`: A1 (tier 1 ≈ 26.8K + memory, tier 2 ≈ 28.7K — the three untouched rules and CLAUDE.md remain the second pass), A3 (vv-principles, coding-elegance; not yet numerical-bug-signatures / algebra-of-record), A4, A5, A6, C5 (settings fix only), D1, D3 (word caps + delegation sentence in the role blocks), E1, H1, H2, H4 (the rule + the page + the role blocks), H5, I1, I2, I6 (as amended), J1–J3, K1, K2, K4 (cores), K7, K8. NOT YET: A7, A8, B1–B4, B5 (the plan floor checklist is still `scratch/_harness_eval/plan-checklist.md` — move it to `docs/development/` in P1b), B6, C1–C4, C6, D2 (session-start.txt carries no Fable-only nudges yet), D4, D5, G1–G7 (acceptance: fresh session), I3, I5, I7, I8, J4–J6, K3, K5, K6.

## Durable lessons from this session (for the evidence page at P1b, not for a rule)

- The harness's agent REGISTRY is fixed at session start; a new agent directory is invisible until a fresh session (a headless `claude -p` session sees it). AGENT.md CONTENT loads fresh per dispatch.
- A sub-agent that emits one Bash call of more than ~10K output tokens is killed by the stream watchdog; chunked writes of ≤ 5 KB per call are the working shape (five agents, zero kills, after two kills the day before).
- `git mv` of a page breaks its relative `:doc:` links (now resolved from the new directory) and every `:doc:` citation of its old name; both are `-W` reds, both mechanical.
- MyST under `-W`: a `###` directly under `#` is an error (`myst.header`); docutils refuses any line over 10 000 characters — a verbatim single-line paragraph needs a content-neutral soft wrap.
- `Write(<path>)` in `settings.json` `permissions.allow` is inert; only `Edit(<path>)` rules match file-editing tools (the harness now says so at startup).
- The three warning classes were found by the build, not by any pre-check — the docs build IS the anchor checker, which is what I6/I3 predicted.

## Resume surface

This section → Part VIII → Part VII. Memory: `project_harness_context_budget.md`. Branch `docs/development-substrate` @ `72006892`. Scratch (untracked): `scratch/_harness_eval/` (drafts, `p1_spec.md`, `p1_gate.sh`, `probe_questions.md`, `p1_commit_msg.txt`, `backup/`, `worklog_2026-09-05.md`).

## T1 and T2 executed (2026-09-20, the post-compaction session)

Five reviews were dispatched in parallel under the new role blocks and skill cores (qa on the rule cores, qa on the skill cores, elegance-enforcer on the generator, a scripted verbatim audit, a qa mutation battery in its own worktree); every return honoured the word cap and the `NEEDS:` block. Reports: `scratch/_harness_eval/review/{qa_fidelity_rules,qa_fidelity_skills,elegance_generator,verbatim_audit,qa_generator_mutation}.md`; the repair scripts are beside them (`t1_*.py`). Everything below was fixed in the SOURCE and regenerated.

| charter step | finding `[M]` | repair |
|---|---|---|
| 1 fidelity to the rulings | roles, invariants, depth 3, parent-review independence, return contract all present; the workflows rule carried two arrow chains (the articulation rule's own tell); the workflows page said "2.1.2xx"; the harness page's cost sum omitted CLAUDE.md | W1/W2 rewritten as phase lists; version made exact; sum corrected |
| 2 distilled cores vs originals | plan-authoring 61 of 67 units in full, 4 lost (distillation obligation; §1 phase template; the inherited-scope rider; §6d's AST edge census), 1 unsafe merge (§6b spelling 12 fused two mechanisms under a check that reaches one), two one-token check defects (`gh issue --search`; `:meth:` dropped), appendix census 110 not reproducible; coding-standards 49 of 49; vv-principles 79 of 80 top level but **8 of 28 test-design sub-items** (Mode 8 at 1 of 11; the only two false-RED guards dropped), `allowed-tools: Bash` lost from the front matter; coding-elegance 98 of 99 (the miss is vacuous) | all four losses restored (the template on the evidence page, linked from MEANS-IS-HYPOTHESIS); (12) split into (12) and (12b); both tokens fixed; appendix rewritten under one stated convention and recounted by script (110 reproduced: §6b 21, §6d 7, §10 9); Mode 8 restored as nine named classes; Modes 9–12 sub-items restored including both false-RED guards; #13 fifth disguise, #17 riders and a new check (i), #29 (g), #35 companion; the DECODER rule into the instrument-doctrine skill X2(7); `allowed-tools: Bash` restored. Budgets: plan-authoring 9000 → 9600 (≈9278 after), vv-principles 9500 → 11000 (≈10706 after) |
| 3 verbatim move | 0 content lost across six pages (plan-authoring 411 of 414, coding-standards 143 of 144, lessons 2024 of 2109, vv-anti-patterns 1304 of 1310, test-design-modes 148 of 163, coding-elegance 253 of 260 lines matched; every unmatched line is text the move ADDED); reverse: 64 of 64 lesson bodies, 94 of 94 surprise rows, 36 of 36 anti-patterns | none needed |
| 4 the generator | `--check` had NO READER (no CI, no test; the Sphinx wiring runs write mode and discards the output; two "in CI" sentences present-tense-false); the slug re-derived MyST's (409 regex anchors vs 396 real: fenced `#` lines minted phantoms); nothing constrained a target to `.claude/`; budget optional in code; crashes past the problems channel; the BEGIN line parsed its own format string; a dead `relative_to` arm; per-link re-parse. Battery on `9818b5ec`: the three claimed arms bite; blind to an agent budget, to a marker-less stale block (the prescribed repair duplicated it and was then blessed), to non-`.md` and same-file links, to orphans, to a missing source (traceback) | generator rewritten: MyST's own `compute_unique_slug` over a markdown-it parse (396 anchors, equal to the build's); `HEADING_ANCHOR_LEVELS` owned by the generator and imported by `conf.py`; manifest entries validated (keys, source exists, target under `.claude/`); budget required on every kind, role blocks included (300 each); marker-less block and END-without-BEGIN are named problems; general link regex; `os.path.relpath`; cached anchors; orphan detection; `tests/test_harness_generated.py` (foundation) is the reader. Battery re-run on the revision: **12 of 12 arms bite**, each with a named message, control green (`scratch/_harness_eval/review/mut2_arm*.log`) |
| 5 the eleven citations | the old `development` page on `main` was the git-workflow page, so every citation already pointed at a page without the section it named; three new targets still lacked the quoted title ("Structural independence applies above the trusted-library line" is a section of the `algebra-of-record` skill; "Why this design" was a section of the retired handoff skill; "Unify after two instances" is Pattern 6) | 5 citations now `:ref:` the corpus doctrine page (`verification-principled-equivalence` ×2, `verification-structural-independence` ×2, `verification-1g-degeneracy`); 3 name coding-elegance Pattern 6; 1 names the algebra-of-record § "The bifurcation pattern"; the structural-independence pair also names the skill's section by its real title |
| 6 frontmatter matrix | matches H5 (Agent on six Key, `omitClaudeMd` on three Support, instrument-doctrine on four); two live pointers to the retired skill in `algebra-of-record`, four present-tense `DISPATCH_REQUEST` instructions in the method-implementer body | all rewritten to the `NEEDS:` contract and the implementer's own `Agent` grant |
| 7 Support-brief paragraphs | — | DEFERRED to T5 by the protocol table (it is the planted task; the measurement needs them absent now) |
| 8 gates | generator `--check` 0 problems / 0 drift; `sphinx -E -W` 0 warnings rc 0 (run 1, before the foundation marker; run 2 on the final tree: `t1_sphinx2.log`); `dead_references` 0 of 66; 422 passed under `-O` across the five harness-facing modules; `git show --stat 72006892 | grep scratch` = 0 | — |
| T2 | `.claude/hooks/log-instructions-loaded.sh` registered under `InstructionsLoaded` (matcher `*`), logs the whole payload one line per event to `scratch/_harness_eval/instructions_loaded.log`; smoke-tested with two events; live from the next session start | — |

Agent `NEEDS:` rulings: the qa's proposed vv-principles anti-pattern "two mechanisms, one check" is not a V&V item (K6); it is the distillation rule, now a bullet in `docs/development/harness.md` § Adding or changing. The 2026-09-04 draft the qa could not obtain is `scratch/_harness_eval/plan-authoring.core.md`; not needed, the appendix is recounted by convention. The L25 row was REMEDIED by extending the lessons index line. The mutation qa's two owed items landed as #17 (i) and the rider on #17 (a).

Durable lessons for the evidence page at P1b: an instrument with no reader (`--check` existed for a day with nothing running it; the docs and the manifest said "in CI" while the repo has no CI); a distilled clause's check must reach every mechanism its text names; `git -C <dir>` resolves a relative path argument inside `<dir>`; the prescribed-repair arm (a gate whose message names a fix owes an arm that applies it and re-checks).

Always-on block after T1 `[M]` chars/3.6: plan-authoring 9.3K, coding-standards 4.9K, instrument-doctrine 1.0K, articulation 0.7K, workflows 1.2K; untouched rules and CLAUDE.md unchanged; ≈27.2K plus MEMORY.md. Session-start batch: vv-principles 10.7K, coding-elegance 8.7K, instrument-doctrine 2.0K, lessons 3.1K, briefing ≈5.5K: ≈30K. Both above the P1 figures by what the four reviews found missing; the T4 targets (memory files ≤ 25K, Key dispatch floor ≤ 75K) stand.

### T1, second pass (2026-09-20, same session; the user ruled "no effort spared")

The first pass could not review two things: the restorations the main agent wrote in response to it, and the pointers from the untouched rules into surfaces the move changed. Four reviews closed that: a qa recount of the restorations against the originals, a qa audit of the lessons-retirement claim with a tree-wide pointer census, an independent qa read of the three new rules, the workflows and harness pages and the nine role blocks against Part VII, and the elegance-enforcer resumed by name on the rewritten generator. Reports: `scratch/_harness_eval/review/{qa_restorations,qa_lessons_retirement,qa_rulings_fidelity}.md` and the "Re-review at 0d376e4f" section of `elegance_generator.md`.

| review | finding `[M]` | repair |
|---|---|---|
| restorations vs originals | 23 of 30 faithful; Mode 11's surrogate clause lost its qualifier and its two-sided mutation check; four modal words dropped; two Mode 8 narrowings; the appendix's prior total (110) was not measured and §10 changed 10 → 9 silently; #17 (i)'s `[M]` named no instrument; the (12b) grep was never run | all restored or re-worded; the appendix carries a `[REMEDIED]` (the pre-review core counted 108 under the stated convention); #17 (i) cites the battery arm; the two new #17 riders gained an evidence entry (`vv-anti-patterns.md` § 2026-09-20 the prescribed repair); L25's body gained a dated addendum; (12b) validated tree-wide: 11 `isinstance` doors on space types in `orpheus/` + `tests/`, the founding case's `SphericalHarmonicSpace` door among them (`tests/numerics/test_frame.py:171`) |
| lessons retirement (20 retired, 44 kept) | 12 of 20 fully carried, 8 partial, 0 lost; 44 of 44 surviving lines faithful, 44 of 44 anchors resolve; two disposition defects (L21's note contradicted its body; L15's ban lost its justification); 9 stale pointer sites at 7 targets, 5 of them new | the ruling's own logic applied: an uncarried mechanism gets a clause in the rule that owns it — the commit-message backtick hazard (L61 b) as a process-discipline section; retire-leaves-first and a move's residues (L20, L34) as coding-standards items 24 and 25; the closure-seed exactness rule (L27) on vv-principles Mode 7; the container census (L55) on plan-authoring VALIDATE-THE-FILTER; the staleness trigger (L8) on nexus-tools; the compute-programmatically check (L4) on CLAUDE.md; L21 restored to the index; L15's line carries the `SweepCoefficientCache` ban; both notes corrected in place; all 9 pointers repointed to evidence anchors; the 50 obsidian `[[lessons-Lnn]]` cross-references on the evidence pages converted to links (MyST rendered them as literal text) |
| rulings fidelity (64 elements) | 56 of 64 encoded; the workflows rule carried the pre-restructure 71K; "2.1.259" was unsupported (the 2026-09-05 probe ran on 2.1.261); the page said W1-P3 is dispatched by the orchestrator where the rule says the parent; W2-P3 had one reviewer in the workflow texts and two in the role blocks; the instrument-doctrine rule claimed the cores cite X1–X4 by ID (1 citation); harness.md's rule shape (check, tell, evidence link) was not met by the three analysis-distilled rules; two metaphors and one undefined term; the K1/K2 caps (800, 600) exceeded (1003, 708) | all fixed at the source: 29.5K `[M]`; 2.1.261; "the parent"; W2-P3 = qa and elegance-enforcer in parallel (a fix is where a patch replaces the structural repair); the ID claim weakened to what is true (adding the IDs to the cores is a P1b item); harness.md's rule shape now says where check/tell apply and that an evidence link comes with the first founding case; "fence" → "restriction", "prophecy" → "claim", "three-leg verdicts" defined; the caps annotated `[M]` in K1/K2 (the manifest budget is the instrument) |
| generator re-review | PASS with nits; all 8 first-pass findings remedied; one concern (the stale-block sentinel keyed on an agent page's H1, which a docs edit could reword, reopening the duplication) and four nits | sentinel keyed on the invariant " — role block" suffix; orphans = any stamped `.md` under `.claude/` at column 0 (one definition); `conf.py` states why the cap lives in the generator; N9: a generated file's link to a manifest source is written to that source's generated copy (the harness copy already in context); battery re-run: markers deleted with the H1 reworded is caught; orphans in a fresh directory and a stray marker are caught |

Also this pass: `git worktree remove` of the leftover `nexus-workspace-wiring` checkout (its branch was merged; the elegance-enforcer reconciled the three memory files it held — one kept, two sections merged, its pre-distillation index dropped); the enforcer's and the two qa agents' memory updates are theirs, committed with this landing.

Gates after the second pass: generator `--check` 0 problems / 0 drift; 422 passed under `-O` across the five harness-facing modules; `sphinx -E -W` and `dead_references` recorded in the landing commit. Sizes `[M]` chars/3.6: plan-authoring 9.4K, coding-standards 5.2K, workflows 1.3K, instrument-doctrine 1.0K, articulation 0.7K; vv-principles 10.9K, coding-elegance 8.7K, instrument-doctrine skill 2.0K, lessons 3.3K. Always-on block ≈27.8K plus memory ≈2.5K; session-start batch ≈30.4K. Budgets: plan-authoring 9600, coding-standards 5300, vv-principles 11000, lessons 3400.

A third, narrow check (`qa_pass2_clauses.md`) read the second pass's own eight clauses and two index lines against their lesson bodies: 5 of 10 faithful at first; the other five were rewritten to the bodies' own words (item 24 defines a leaf as zero production callers and keeps the audit table at `.claude/plans/<plan>_dependency_audit.md`; VALIDATE-THE-FILTER carries L55's separator-anchored check; the staleness trigger is the body's two symptoms after a major move; Mode 11 keeps its "only catchers" scope and "structurally blind"; L15's ban is on the NAME), and `git_workflow.rst`'s published `git commit -m` recipe became the `-F` form.

P1b items this pass adds: cite X1–X4 by ID from the domain cores; the 21 glyphs in the three untouched rules (K4 says retired) go when those rules come under generation; 46 `[[lessons-Lnn]]` wiki-links in `.claude/agent-memory/**` and `.claude/plans/**` are agent- and plan-owned (20 point at retired numbers) and are re-pointed by their owners.

## T4 executed (2026-09-20, the session restarted on the branch)

Claude Code 2.1.278, restarted on `docs/development-substrate` @ `6e628230` (the `InstructionsLoaded` log carries that restart's `session_start` events at 11:27:04Z), then `/compact` at 11:30:23Z. Three instruments, each with both legs stated. Probe transcripts: `~/.claude/projects/-Users-rodrigo-git-nuclear-ORPHEUS/{c03a470c…,bb6438d4…}/subagents/` (the F10 fixture, old and new substrate) and `54d619d4…/subagents/` (the two dispatches from this session); the headless parent's output is `scratch/_harness_eval/t4_headless_probe.out`.

**1. The always-on block** (CLAUDE.md, the 8 unscoped rules, MEMORY.md — the `instructions` attachment of a Key-shaped dispatch).

| instrument | old substrate (`main`) | new substrate (branch @ `6e628230`) |
|---|---|---|
| keep − omit, first-turn input tokens (`cache_creation + cache_read + input`) on the F10 fixture: throwaway `_probe-keep` / `_probe-omit`, haiku, one-line system prompt, 0 tools, dispatched by a headless `claude -p --model sonnet`; `[M]` sessions c03a470c and bb6438d4, both 2026-09-20 | 107 485 − 36 519 = **70 966** over 249 137 chars in 7 files (CLAUDE.md 9.6K, plan-authoring 159.8K, coding-standards 36.7K, process-discipline 12.0K, nexus-tools 9.3K, delegation 4.5K, MEMORY.md 17.2K) | 74 822 − 44 167 = **30 655** over 108 284 chars in 10 files (CLAUDE.md 10.1K, plan-authoring 33.8K, coding-standards 18.6K, process-discipline 12.7K, nexus-tools 9.5K, delegation 4.6K, workflows 4.4K, instrument-doctrine 3.4K, articulation 2.3K, MEMORY.md 8.9K) |
| chars / 3.6, the generator's estimate | 69.2K (64.4K without MEMORY.md — the "≈65K" the target was set from) | 30.1K |
| `/context`, the "Memory files" category | never recorded on the old substrate | 41.8K, read after `/compact` in this session |

Target ≤ 25K: **missed by 5.7K on the API tokenizer** (−57 % from 70 966). The overrun sits entirely in the four files P1 left untouched by ruling (K3, P1b): CLAUDE.md + process-discipline + nexus-tools + delegation = 36 992 chars ≈ 10.5K tokens; the five generated cores are 62 431 chars ≈ 17.7K and MEMORY.md is 8 861 chars ≈ 2.5K. Three facts about the instruments: (a) the API tokenizer reads both substrates at 3.51–3.53 chars/token (`[M]` 249 137 / 70 966 and 108 284 / 30 655), so `CHARS_PER_TOKEN = 3.6` in the generator is calibrated and its token budgets mean what they say; (b) `/context`'s category figure over-reads the same block by 36 % (41.8K against 30.7K), so G1's instrument is the keep − omit probe, not `/context`; (c) the harness strips the `<!-- GENERATED … -->` stamp before injection (`[M]` workflows.md is 4 564 chars on disk and 4 350 in the attachment; the same ≈215-char delta on every generated rule), so the stamp costs no context.

**2. The dispatch floor.**

| probe | old | new |
|---|---|---|
| Key-shaped, the F10 fixture (`_probe-keep`: haiku, 0 tools, headless parent), first-turn input | 107 485 | **74 822** — under the 75K target. The harness-fixed block itself grew between the two runs: the `omit` arm reads 36 519 then and 44 167 now; the four claude.ai connector servers `/context` lists (14 MCP tools) account for ≈3.3K of the 7.6K, the rest is unattributed `[R]`. On the F10-era fixed block the new figure would read ≈67.2K `[R]` (74 822 − 7 648, assuming the blocks add) |
| Key-shaped, the G2 fixture (`general-purpose`, haiku, 0 tools, dispatched from the interactive session), `subagent_tokens` as the task notification reports it | 113 006 | 88 679; first-turn input 82 746. Not the same instrument twice: the notification figure now includes the `SubagentHandback` turns (`tool_uses` 1 against 0 then), and `general-purpose` is dispatched with every tool schema (`Tools: *`) where the F10 probe is dispatched with none — that, not the substrate, is the 8K between 74 822 and 82 746 |
| Support: `_probe-omit` (F10 fixture) and `explorer` (`omitClaudeMd`, two preloaded skills, haiku, from the interactive session), first-turn input | 36 519 / not measured | 44 167 / 37 717. The explorer's transcript carries no `instructions` attachment (`omitClaudeMd` holds); its two preloaded skills arrive as user turns, 11 398 chars |

Answers: both `keep` arms and the `general-purpose` probe answered YES to the CLAUDE.md heading and to "ugrep 7.5.0", NO to the two retired headings (plan-authoring's "Surprise log — each clause traces to what produced it"; coding-standards' "Retire as you go (aggressive retirement)"), NO to "Four index disciplines" (the new index says "Index disciplines") and quoted the new Active-work bullet: the dispatched context is the branch's substrate, not a stale snapshot. Every `omit` arm answered NO to items 1–6 and YES only to the harness-injected `Co-Authored-By` line.

**3. The `InstructionsLoaded` log** (`scratch/_harness_eval/instructions_loaded.log`; 27 events in three sessions). `session_start` (03a2f2b6, the restart; bb6438d4, the headless parent) and `compact` (54d619d4) each log exactly 9 files — CLAUDE.md and the 8 unscoped rules (workflows, delegation, instrument-doctrine, plan-authoring, nexus-tools, process-discipline, coding-standards, articulation), `memory_type` `Project` on every event. Not logged: MEMORY.md (the auto-memory index is not an `InstructionsLoaded` event); the path-scoped `vv-testing.md` (no `tests/**` path was touched); and sub-agent dispatches — 0 events after four of them (two Key-shaped with the `instructions` attachment present, two Support with `omitClaudeMd`), which `[REFUTED 2026-09-20]` the T2 row's "one dispatch per role" expectation while its main-agent half holds. No retired file appears: `git diff --diff-filter=D --name-only main HEAD` names one file, `.claude/skills/subagent-handoff-protocol/SKILL.md`, a skill and never a rule. Hook artefact: concurrent events interleave inside a line (the per-event `printf` + `tr` is not atomic across processes), so the log is parsed by a regex on `"file_path":…"load_reason"`, never by line.

**Verdict.** The Key-shaped floor target is met on the fixture that set it; the memory-files target is missed by the files the ruling deferred to P1b, not by the cores; the hook measures the main agent and nothing else. Nothing here changes the order T5 → T6 → T7. The throwaway agents were removed (`git status --porcelain .claude/agents/` is empty).

## T5 executed (2026-09-20, the same restarted session)

The planted task: the "Support briefs" paragraph in the six Key role blocks, written by the main agent with one Support dispatch (explorer, a census brief written to the template with a PLANTED positive control) and one Key dispatch (qa, W4 claim verification, resumed by name on the repaired tree). Reports: `scratch/_harness_eval/t5/{explorer_brief_census,qa_support_briefs}.md`; Sphinx logs beside them.

**The Support dispatch.** Brief: W4-P0, the template's five parts, four rules under "Rules that apply to you" (ugrep silent zero; a NAMED control, `delegation.md` § "Briefing a literature pull"; predicate and tree per count; no Nexus, and NEEDS: on a missing tool), return contract 300 words + file + NEEDS:. Return `[M]`: 321 words (7 % over the cap); the control validated on the first run and its blind spot named (two lines of the control section matched no net, found by reading); "47 entries in 30 of 94 files under <roots>" with the predicate; 8 of 8 hidden members found (CLAUDE.md:101, the template, invariant 4, the explorer role block, nexus-tools:113, the behavioral-auto-regression skill, L28, A-BRIEF'S-METHOD) plus three the key had not listed (L38, L50, L12); the tree split HEAD / working tree when it changed under the census; Nexus untouched (19 tool uses, Bash/Read/Write). Two findings of its own: L12's body cited a template section that no longer exists (repaired in place, `[REMEDIED]`); the shell's ugrep is 7.8.4, not the 7.5.0 the rule measured (`[M]` the fixture reproduces the silent zero on 7.8.4; noted in the rule).

**The Key dispatch.** qa returned 12 findings under the cap with a NEEDS: block; it dispatched its own zero-tool explorer probe to test the paragraphs' premise. `[M]` F1: `omitClaudeMd` drops CLAUDE.md, the rules and the PROJECT memory index, not the agent's own memory index (`.claude/agent-memory/explorer/MEMORY.md`, 14 311 B) or its preloaded skills — six role blocks and the template had said "no memory index" from F10's number, which measured the project index only. F3/F7/F8: the paragraphs restated four template items and dropped the two "Always" ones. F4: CLAUDE.md:101 drifted ("two or three rules"). F5: invariant 4's "sees no project rule" was false for the haiku categorisers the roles table calls Support. F6: cross-domain-attacker named nowhere. F9: delegation.md's literature section is fuller than the template's clause, no cross-reference. F10: the manifest comment said what a budget counts, not how it is set. F11/F12 (facts, no repair): the budget was the only objecting instrument and the raise removes its objection — it is a size gate, the review is the quality gate; the F10 probe fixture carries no role block, so the T4 floor cannot see a role-block change (≈0.2K per Key dispatch here). Repairs, all in `docs/development/` and regenerated: the template's "Rules that apply to you" line is the ONE definition (three `omitClaudeMd` agents; paste the line in; Always: L28, read-only or not, `python -O -m pytest`, L38, L50; census: the ugrep remedy, a control per shape, exclusions; literature: delegation.md § in full; design review; Nexus fallback; pointers verbatim) and the return contract carries L12; the six paragraphs point at it and paste it, naming only the founding exposure; invariant 4 and the roles row name the three `omitClaudeMd` agents; CLAUDE.md:101 points at the template; the manifest states how a budget is set; the harness page gains "a claim about what a dispatch receives is measured by dispatching". Rulings on the qa's NEEDS: the wording (adopted as above); its vv-principles drop-in is a harness lesson, not a V&V item (as K6); it may write its own agent memory. Confirmation on the repaired tree: 11 of 12 closed, one half, with six residues R1–R6 (`[M]` R1: the manifest comment's "next hundred above the measured size" was false for 4 of 8 text budgets and `--check` asserted only `size <= budget`). Residues closed by the orchestrator without a third qa round: the convention restated as "a round figure at most SLACK_MAX = 400 tokens above the measured size" and ASSERTED — `generate_harness.py` now fails a budget with more slack (`[M]` articulation's budget set to 8000 on a copy-aside manifest reddens `--check` with the SLACK line; restored byte-identical); method-implementer's May-call names cross-domain-attacker; the census clause adds the Python re-run of a completeness claim; the three Support role blocks point at the template instead of restating it; the lessons index and the template agree on "verbatim, in a code fence". The qa wrote its own lessons to `.claude/agent-memory/qa/` under ruling.

**Surprises by root** (2, both REPEATS of clauses that exist): (1) root L38 — the orchestrator edited the censused tree while the census ran ("about to be added" in the brief; added five minutes later); the explorer split its count by tree; the writer's obligation is now the template's L38 line. (2) root §2 [M]-SCOPE — "no memory index" written from a measurement of a different quantity; logged as the surprise-log entry `2026-09-20 omitClaudeMd memory index`, linked from the clause.

**Contracts.** Word caps: explorer 321/300, qa under 400 and under 250 on resume; both returns ended with NEEDS:; qa's NEEDS were acted on and it was resumed by name (invariant 3). Budgets raised in the same commit as the text: the workflows rule 1300 → 1400 (`[M]` 1307 before the repair), the six Key agents 300 → 500 (`[M]` 386–449). Gates: `generate_harness --check` 0/0; `test_harness_generated` passes under `-O`; `sphinx -E -W` run 2 (the repaired tree before the residues) rc 0; run 3 (the committed tree) rc 0 (`t5/sphinx.log`); Nexus `dead_references` on the rebuilt graph 0 dead of 66 checked.

## T6 executed (2026-09-20, headless, both substrates)

The recall test K5 ruled: one past qa review whose findings a landed fix confirmed, re-run from the same brief on both substrates, scored against the known finding set with a denominator. Fixture, scoreboard and every report: `scratch/_harness_eval/t6/` (`brief.md`, `rubric.md`, `scoreboard.md`, `<arm>_findings.md`, `<arm>_report.txt`).

**Fixture.** The review round of CS4c step 5: `scratch/cs4c_step5_qa_report.md` reviewed the test re-keys of `8d432a5d` (37 files; 24 test files, +2569/−446) and its ten findings F-1 to F-10 were confirmed by `caeb995d` (nine landed as test changes; F-3 became the vv-principles #22 widening). Four detached worktrees at `8d432a5d` under `ORPHEUS_t6/` (a sibling of the repo, so the main checkout and its auto-memory are untouched), each with one substrate overlaid — CLAUDE.md and `.claude/` from `main` @ `7e13621b` or from the branch @ `df3e0f31` — crossed with two reviewers: the Opus-pinned `qa` agent dispatched by a headless sonnet parent, and a Fable main-session review (`claude -p --model fable`). The SAME brief file in all four: W1-P3, the diff as a patch file, the implementer's summary and the verification plan (both still in `scratch/`), read-only, in-process mutations only, the history beyond the two commits forbidden, Nexus absent by construction (`--mcp-config {} --strict-mcp-config`), and a pytest wrapper that removes the venv's editable finder so each worktree's own `orpheus/` is under test (`[M]` `orpheus.__file__` per worktree; 82 tests collected). `--dangerously-skip-permissions` on all four so the two substrates' allow-lists are not a confound. The qa arms needed `CLAUDE_CODE_PRINT_BG_WAIT_CEILING_MS=0`: print mode terminates background sub-agents at 600 s (`[M]` both first attempts cut there; relaunched from a reset scratch folder). Contamination checked: the qa's agent memory on both substrates carries no step-5 entry; no arm ran a history command (`[M]` 0 of 4 transcripts; the `#452` both old arms cite is in the old and the new vv-principles skill, item 36).

**Scoreboard** (a finding is reproduced when the arm names the same test or site and the same mechanism; scored by reading; ½ = half the mechanism):

| known finding | old_qa | new_qa | old_fable | new_fable |
|---|---|---|---|---|
| F-1 reciprocity row pairs a signed pairing with a squared norm | ✓ (cites the round from the rule's row) | ✓ | ✓ | ✓ |
| F-2 (n,2n) sibling row compares 0.0 with 0.0 | ✓ | ✓ | ✓ | ✓ |
| F-3 conjugation identity shares its middle factor | ✗ | ✓ | ✓ | ✗ |
| F-4 displaced guard arm lost its only witness | ✗ | ✗ | ✗ | ✗ |
| F-5 C6 static half runs nowhere, narrowed | ✓ | ✓ | ✓ | ½ |
| F-6 coverage shift: k∞ no longer witnesses `apply` | ✓ | ✓ | ✓ | ✓ |
| F-7 `hasattr`-guarded conjunction | ✓ | ✓ | ✗ | ✓ |
| F-8 G5.1 ships 3 of 4 clauses | ✗ | ✗ | ✗ | ✓ |
| F-9 `_PACKAGE` cwd-relative | ✓ | ✓ | ✓ | ✗ |
| F-10 plain/lift row is a declared wrapper comparison | ✗ | ½ | ✓ | ✓ |
| **k of 10** | **6** | **7.5** | **7** | **6.5** |
| **k of 9, F-1 excluded** | **5** | **6.5** | **6** | **5.5** |
| findings beyond the known set | 10 of 16 | 7 of 15 | 9 of 16 | 6 of 13 |
| first-turn input tokens | 233 432 (qa); 134 703 (parent) | 111 678 (qa); 86 140 (parent) | 136 542 | 87 979 |
| assistant turns / tool uses | 123 / 69 | 121 / 73 | 81 / — | 112 / 76 |
| report words (cap 600) | 503 | 693 | ≈620 | 654 |

**Reading.** F-1 is CONTAMINATED on every arm and is excluded from the denominator: the always-on plan-authoring rule carries its exact mechanism on both substrates (old: the 2026-09-05 surprise row, which names the round, the report path and the remedy — old_qa cited all three; new: the INERT-WITNESS clause). **k of 9 with F-1 excluded: old_qa 5, new_qa 6.5, old_fable 6, new_fable 5.5.** F-4 (the displaced guard arm) was found by no fresh arm and F-8 by one; both are findings the original reviewer reached from the plan's own text (the arm's message, the fourth clause), the part of the fixture a fresh reviewer reads least. 3 of 4 fresh arms (new_qa, old_fable, new_fable) flag one demotion the original review missed: `TestBitIdenticalExtractionP0::test_isotropic_kernel_bit_identical_to_legacy_verbs` replaced the solver's cached pair by a test-minted pair while its comment credits the production pair (new_qa proved it with a descriptor read counter: zero production reads). Verify on `main` before filing; the class has since been reworked. One run per cell is one draw (K5's own caveat), so the substrate comparison is bounded, not decided: the cores lose no recall that this fixture can see, and the Fable main-session review costs 48.5K fewer first-turn tokens on the new substrate (`[M]` 136 542 → 87 979).

**Contracts.** Every arm ended with `NEEDS:`; three of four arms exceeded the 600-word cap (≈620, 654, 693 words); old_qa kept it (503); no arm edited a tracked file (`[M]` `git status` per worktree). Worktrees removed after the archive; `git worktree prune`.

## T7 executed (2026-09-20)

Filed before the merge so this record carries the numbers: J4 → sphinxcontrib-nexus#94 (eager MCP tool schemas against the deferred hint, with the per-tool sizes); J6 → sphinxcontrib-nexus#95 (the generic ontology-driven directive, VII.4); J5 → ORPHEUS #477 (the umbrella queue, `module:docs`, carrying the proposal IDs and the P1b list); J1 → the landing comment on #308. Merge: `git checkout main && git merge --ff-only docs/development-substrate` — this commit is the tip, so the check is `git merge-base --is-ancestor 87f14835 main`; then `git push origin main` and `git branch -d docs/development-substrate` (the branch was never pushed, nothing to delete remotely). CI: this repository has none (`[M]` `ls .github/workflows` prints nothing; `gh run list` prints nothing), so the row's "read CI after the push" instrument does not exist here; the gates that stand in for it ran on this tree — `generate_harness --check` 0 problems / 0 drift, `test_harness_generated` under `-O`, `sphinx -E -W` rc 0 (T5 run 3), `dead_references` 0 of 66. The campaign memory's terminal state becomes "merged" at merge time (process-discipline); the queue is #477.

## Transition and evaluation protocol (ruled 2026-09-20; execute in this order)

The harness's asymmetry is the method: rules and CLAUDE.md are a SESSION-START snapshot, skills and AGENT.md load fresh per invocation or dispatch, and rules are read from the WORKING TREE, not from `main`. So the transition is staged, and the new substrate is evaluated on the branch before anything is merged.

| step | session | what | instrument, and what it must read |
|---|---|---|---|
| T1 | this one, after compaction | run the review charter above. Its dispatches already run on the NEW skill cores and role blocks (live per dispatch) — that is half the change under real load, pre-merge | qa / elegance-enforcer / general-purpose returns honour the return contract (word cap, `NEEDS:`); any finding fixed in the SOURCE, regenerated, gates re-run |
| T2 | this one | add the `InstructionsLoaded` hook to `.claude/settings.json` (logs every instruction file the harness injects, and why — `code.claude.com/docs/en/hooks#instructionsloaded`), writing to `scratch/_harness_eval/instructions_loaded.log` | the first new-substrate session logs itself: main agent + one dispatch per role (Key with rules; Support with `omitClaudeMd`) `[REFUTED 2026-09-20]` for the dispatch half — 0 events after four dispatches; the main-agent half holds (9 files per `session_start` and per `compact`); see "T4 executed" |
| T3 | — | compact; RESTART ON THE BRANCH `docs/development-substrate` (unmerged) | the new rules load from the working tree; `main` stays the rollback |
| T4 | fresh, on the branch | `/context` → memory-files block; the haiku probe (`scratch/_harness_eval/probe_questions.md`, zero tools) → first-turn usage from its transcript; read the hook log | targets: memory files ≤ 25K (was ≈65K); dispatch floor ≤ 75K for a Key-shaped agent (was 113K), ≈36.5K for Support; the log names exactly the generated files and no retired one. **`[M]` 2026-09-20 executed — "T4 executed" below: block 30 655 (missed, by the untouched files); Key floor 74 822 (met); Support 44 167; the log names the 9 always-on files and no retired one** |
| T5 | fresh, on the branch | one SMALL real task under the new rules — P1b's support-brief paragraphs in the six Key role blocks (ten lines each) — done by the main agent with one Key dispatch and one Support dispatch | contracts honoured; the explorer dispatch gets a census brief with a PLANTED positive control (a known member the filter must find) — does the Support brief carry enough rule for it to validate its filter? Surprises logged by root. **`[M]` 2026-09-20 executed — "T5 executed" below: the explorer validated the planted control first and found 8 of 8 hidden members; the qa's probe refuted "no memory index"; two surprises by root** |
| T6 | headless, both arms | the recall test (K5): one past qa review whose findings a landed fix confirmed (candidate: the #426 review round, `f52877db`; or any `scratch/*qa_report*.md` with a landed fix); check out the pre-fix commit's diff; run the SAME review brief from two `claude -p` sessions — one launched in a `git worktree` at `main` (old substrate; L22: build the graph there or brief the agent to grep), one on the branch — and score each against the known finding set | `k of N` findings reproduced per arm, on both models where the agent is Opus-pinned; a denominator, not an impression. Report both numbers with the fixture. **`[M]` 2026-09-20 executed — "T6 executed" below: of 9 (F-1 excluded, contaminated) old_qa 5, new_qa 6.5, old_fable 6, new_fable 5.5; first-turn tokens old → new: qa 233 432 → 111 678, Fable main session 136 542 → 87 979** |
| T7 | fresh, on the branch | merge: `git checkout main && git merge --ff-only docs/development-substrate`; push; delete the branch; comment #308 with the hash; file J5 (umbrella) + J4/J6 (nexus) | `git merge-base --is-ancestor` true; CI read after the push (process-discipline: baseline a red CI before adding to it). **`[M]` 2026-09-20 executed — "T7 executed" below: sphinxcontrib-nexus#94, #95; ORPHEUS #477; #308 commented; no CI exists in this repository** |
| T8 | fresh, on `main` | resume the campaign queue (P1b → skills split → tools/hooks → root layer); the surprise log tagged by root on the first two campaigns is the ongoing detector | surprises per root; a recurring root is a TOOL, not a paragraph |

What NOT to do: judge adherence from one session impressionistically; run full-vs-core on one prompt once per arm (indistinguishable from run-to-run variance — the reason K5 replaced the A/B); merge before T4–T6 have numbers.

Rollback at any step: `git checkout main` (old substrate loads at the next restart); the old MEMORY.md is `scratch/_harness_eval/backup/MEMORY.md.20260920`.
