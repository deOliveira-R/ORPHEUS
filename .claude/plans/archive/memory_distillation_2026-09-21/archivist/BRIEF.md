# Brief — memory distillation, owner pass: archivist

**Workflow:** memory distillation (the harness campaign, `.claude/plans/harness_context_budget.md`, ⏸ COMPACTION POINT #11). **Phase:** PROPOSE. You are the OWNER of `.claude/agent-memory/archivist/`; you propose, the orchestrator evaluates and applies, one commit per owner. **You edit nothing under `.claude/agent-memory/`, `.claude/agents/`, `.claude/rules/` or `.claude/skills/`: every output is a NEW file under `scratch/_memory_distillation/archivist/`.**

## The ruling (the user, 2026-09-21, verbatim)

> "We will do the memory distillation next for all agents that got bloated memories and any duplicate entries (such as the memory you have to neglect effort and files touched when deciding if something should not be done, which is now a rule. If the rule is not sticking, we need to improve the rule, not by having duplicates)."

Read as three laws: (1) the distillation covers every agent memory that bloated; (2) a memory that restates a rule is a duplicate and retires at the memory, the rule cited in its place; (3) a rule that is not sticking is improved at the rule (its check, its tell, its tier), never shored up by a second copy in a memory.

## The standard (the main memory `feedback_memory_distillation_standard.md`, pasted whole)


The test for whether a memory line earns its place: **"what mistake did I
make, if any, and what did I learn in that session that improved my
behaviour?"** A line that answers it (a failure + the behavioral correction)
stays. A line that is campaign narrative — "in Wave O this and this was done",
"in Phase C bla bla" — is **noise**: *what is Wave O? what is Phase C?* The
name is a pointer into git archaeology that brings nothing forward. The point
was never the campaign; it was the lesson.

**The objective is compression ONLY if it can be done losslessly** (cut prose,
keep the behavioral content + a pointer to the topic file that holds the
detail). **For lossy cases, the objective is DISTILLATION** — find the
meta-lesson that drives multiple failure points and encode THAT, rather than
truncating one instance and losing the rest.

**Meta-lessons and lessons are uplift candidates to `.claude/rules/`.** The
governing question for the whole substrate is *"how to properly behave in this
project?"* A behavioral correction that generalizes across agents/modules
belongs in the always-on floor (`.claude/rules/`), not buried in one agent's
campaign memory. [[harmonize-substrate-by-dimension]] is the sibling rule for
*reconciling* trees; this one is for *distilling* a single tree.

**Process (user, 2026-06-21):** (1) commit all memory first — a checkpoint to
roll back to. (2) Fleet-wide but **one owner at a time** — dispatch the owner
(or a curator) to PROPOSE its distilled memory; evaluate the proposal
carefully before applying; commit per-owner (per-owner rollback points).
(3) Carry meta-lessons forward as rules-uplift candidates surfaced to the user.

**Why (user):** memory exists to sharpen future behaviour. Prose that requires
"deep digging into git history to find out" what a session-internal codename
meant is pure attention tax with zero behavioral payload — and "nothing is
more 'environment' than your context and attention." A 160 KB index loaded
every dispatch is the failure mode this standard prevents.

**CRITICAL — the win is the INDEX slim, NOT file deletion (2026-06-22 lesson):**
The per-dispatch context tax comes from the `MEMORY.md` INDEX (loaded in full
every dispatch). Slimming the index to behavioral hooks banks the ENTIRE win;
the topic `*_closeout.md` files are loaded ON-DEMAND and cost nothing per
dispatch as long as they are not in the index. **Agent-memory files are NODES
in the unified Nexus knowledge graph** — Nexus indexes `.claude/agent-memory`
ALONGSIDE `docs/` and `.claude/plans/`, so a closeout is routinely referenced
from theory pages and campaign plans. Therefore RETIRING a memory file needs
the FULL retirement-audit blast-radius (the graduated `coding-standards.md`
rule: grep `docs/` + `plans/`, which `-W` does NOT warn on) — NOT just the
within-agent-memory wikilink check. Lesson learned the hard way: a finale
`git rm` of 333 closeouts broke 98+ in-repo references (incl. an active
pyright keystone) and had to be reverted. **Default: slim the index (hot
view) + keep the closeouts as cold archive (graph-referenced, un-loaded).**
Delete a memory file only after its full blast radius (docs + plans + graph)
is clean.

**CRITICAL SHARPENING — "slim the index" banks the win only if the index is the
ONLY hot surface (measured 2026-08-03).** The 2026-06-22 rule above is right that
cold topic files cost nothing per dispatch. But it silently assumes the index is
the sole always-loaded file. **The hot surface is actually the TRANSITIVE CLOSURE
of what a loaded file tells you to read.** The moment an index says
*"lessons.md — read FIRST each dispatch"*, that file is hot too, and slimming the
index merely DISPLACES the tax onto it rather than removing it. Measured across the
fleet 14 months later: indices were all healthy (≤200 lines, as designed), while
`archivist/lessons.md` had reached 3820 lines / 254 KB, `archivist` 2808 /
192 KB, `qa` 2572 / 181 KB. Two independent failures follow, and BOTH are silent:
- **Truncation eats the NEWEST.** A default `Read` returns 2000 lines, so an
  append-ordered lessons file past that length drops its most recent — the most
  expensively-learned — entries first. 37 lessons across three agents were
  unreachable, including two written the same day they were measured.
- **Or, if fully paged, a 45–64 K-token bill before the agent does any work.**
There is no good horn. **The fix is structural, not compression:** apply the SAME
hot/cold split one level down — `lessons.md` becomes a ≤400-line digest of
behavioral rules, `lessons_archive.md` keeps the war stories cold. Make the archive
a **pure `git mv`** so the split is lossless by construction and the reviewer can
prove it with one `diff` against `HEAD`.
**Generalised rule: when auditing memory cost, follow every "read FIRST" pointer and
measure the file it names — never stop at the index.**

**During distillation, OFFER each agent the AGENT.md-promotion option (added
2026-06-22).** A distilled lesson is not always memory material: if it is a
STANDING identity-level operating principle — applied on essentially every
task, core to HOW the agent works (e.g. a definition-of-done bar, a standing
review/diagnostic discipline) — it belongs in the agent's always-loaded
**AGENT.md DEFINITION**, not in recalled `lessons.md`. So the distillation
brief MUST offer: *"is any lesson identity-level? If so promote it into
AGENT.md as a standing directive (no war-story, no codenames), leaving the
specific instance in lessons.md with a `→ now in AGENT.md` pointer."* Be
conservative — a thin AGENT.md of sharp identity principles beats a bloated
one; a specific failure→correction stays a lesson. This was MISSED in the
first fleet distillation and run as a separate follow-up pass — bake it into
the distillation itself going forward. Timing bonus: an AGENT.md edit is live
on the agent's **next dispatch** (loaded fresh), whereas rules/MEMORY.md edits
are session-start snapshots that only go live next session
([[harness-context-snapshot-timing]]).

## What is measured (`[M]` 2026-09-21, `wc -l`)

- `lessons.md` 3 256 lines (the hot digest, read every dispatch); `lessons_archive.md` 13 714 lines (cold, sections `L-0NN`); `MEMORY.md` 159 lines (loaded whole every dispatch); 83 topic files (`feedback_*.md`, `425_*.md`, campaign memos).
- References from `docs/`, `.claude/plans/` into `.claude/agent-memory/archivist/lessons.md`: 0. References from `docs/` and `.claude/plans/` into ANY archivist memory file: 0, so the topic files carry no external blast radius either; still, this pass proposes their retirement in the table only and does not perform it.
- A default `Read` returns 2 000 lines: read `lessons.md` in pages (`offset`/`limit`) until you have seen its last line, or the newest entries are the ones you never read.

## Scope of this pass

The two hot surfaces: `lessons.md` and `MEMORY.md`. The archive stays as it is (cold, append-only; you may propose ADDITIONS to it, never edits). Topic files (`*_verification.md`, `active_campaigns.md`, `feedback_*.md`) are not retired here; list the ones you judge archaeology in the table, with the reason, and the orchestrator runs their blast-radius audit separately.

## What counts as a duplicate (law 2)

You hold, every dispatch, the skills `retirement-audit`, `instrument-doctrine`, `nexus-verification`, `nexus-exploring`, `vv-principles`, `algebra-of-record`, and the always-on rules `cardinal`, `articulation`, `instrument-doctrine`, `process-discipline`, `workflows`, `code-search`, `nexus-tools`; touching `tests/**` loads `vv-testing` and `coding-standards`, touching `orpheus/**` loads `coding-standards`, touching `.claude/plans/**` loads `plan-authoring`. A digest entry that restates a check, a tell, an anti-pattern, a mode or a clause of any of these is a duplicate: it retires, and at most one line remains that cites the clause. Adjudicate by READING both sides: open the rule file (`.claude/rules/<name>.md`) or the skill (`.claude/skills/<name>/SKILL.md`) and name the heading or item number in the table; a verdict from memory of the rule is not a verdict. A founding case the rule does not carry is not a duplicate: it is either a lesson that stays (a specific failure → correction) or an uplift candidate (law 3: the rule gains the check or the tell that would have caught it; you propose the clause, the orchestrator edits the rule's source).

## The test for every line

*What mistake did I make, and what correction changed how I work?* A line that answers it stays, cut to its imperative and its `→ LNN` pointer. Campaign narrative (a phase name, a step number, a tracker id, a "⏹ DELIVERED" status, a scratch path that no longer exists) is archaeology: it moves to the archive section that already holds the war story, or to a NEW archive section you write in `archive_additions.md` when no section holds it. A gap that has since closed is STALE: name the landing (the digest's own rule, `MEMORY.md` §1). Several entries that are instances of one meta-lesson DISTIL into that meta-lesson, each instance kept as a pointer only. The standard's figure for the digest is ≤ 400 lines; the count is a reading, never the criterion: distil by the test, then report the count.

## The AGENT.md promotion offer (the standard, 2026-06-22)

Is any lesson identity-level, applied on essentially every task, core to HOW you work (a definition-of-done bar, a standing review discipline)? Propose it as a standing directive sentence for `AGENT.md` in `agent_md.md`, no war story, no codenames, and leave the specific instance in the digest with a `→ now in AGENT.md` pointer. Be conservative: a thin `AGENT.md` of sharp identity principles beats a bloated one. Read the current `AGENT.md` first: its standing directives may already carry what you would promote.

## Deliverables (files; the reply is ≤ 300 words)

1. `lessons.proposed.md` — the distilled digest, complete, ready to replace `lessons.md` verbatim.
2. `MEMORY.proposed.md` — the index with §2 slimmed under the index disciplines (a landed or delivered campaign is ONE line: name, terminal status, pointer; a finding inline there is either a lesson, so it is in the digest, or an open ruling, so it is in a GitHub issue or the archive section; a hook is ≤ 15 words).
3. `table.md` — the audit trail, one row per digest entry (or contiguous line range), columns `lines | entry (first words) | verdict | where its correction lives now`; verdicts: KEPT · CUT (lossless) · DISTILLED into <meta-lesson> · DUPLICATE of <rule/skill § or item> · ARCHAEOLOGY → archive §LNN or `archive_additions.md` §new · STALE (landing named). Then the same for `MEMORY.md` §2. Then the counts: lines before/after per file, rows per verdict (X2).
4. `archive_additions.md` — new archive sections, if any (append-ready, numbered after L86).
5. `uplift.md` — meta-lessons that generalise beyond this agent, each with the rule it would join and the clause it proposes (a `check:` and a `tell:` in the rule's own form), and the founding case it would cite. Empty is an answer.
6. `agent_md.md` — the promotion candidates, or "none".
7. In the reply: the counts, the file list, and a `NEEDS:` block for anything you could not resolve (a lesson whose home you cannot find; a landing you cannot verify against the tree; a rule clause you believe exists but could not locate).

## What NOT to do

Edit or delete any file outside `scratch/_memory_distillation/archivist/`. Distil by character count. Retire a line whose correction no rule and no remaining line carries (the table must show, for every retired line, where the correction lives now). Write a rule's text into the proposal (cite it). Read the archive whole (open a section when a pointer needs checking). Rely on memory of a rule for a DUPLICATE verdict. Use Nexus for this (it is not a graph question; `grep` via Bash is the instrument if you need to check a reference).
