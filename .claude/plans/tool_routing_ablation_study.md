# Tool-Routing Ablation Study — design v0.1

**Status:** DESIGN (not yet executed). Author: main session, 2026-06-14.
**Purpose:** Empirically determine the *minimal* project-instruction dose at which
Claude (Opus 4.8 / Sonnet 4.6 / Haiku 4.5) routes code-exploration tools correctly —
Nexus for structural questions, `grep`/`Read` for known-target point lookups — and to
detect the two failure modes that bracket correct behavior. The result sets the
*dose* of the future `.claude/rules/nexus-tools.md` instead of guessing it.

## 0. Background / why this study exists

- Confirmed 2026-06-14: the standalone `Grep`/`Glob` tools are GONE for all three
  models; Bash discourages only `cat/head/tail/sed/awk/echo` (not `grep/find/rg`);
  no "always Grep" mandate survives. All three report freedom of choice.
- Therefore the OLD bias (always-grep, directive-driven) is dead. But **training
  priors** remain and can bias behavior with no directive at all. Two priors matter:
  - pretraining corpus → reflexive `grep/cat/find` ("under-Nexus").
  - agentic RLHF that rewarded "use the dedicated/special tool / delegate" →
    a transferable "reach for the fancy tool" habit that can land on Nexus as
    theater ("over-Nexus").

## 1. Target behavior (the gold standard)

A response is correctly routed when the FIRST substantive exploration call matches
the question's native structure:

| Question shape | Correct first tool | Wrong (failure mode) |
|---|---|---|
| Callers / dependents / blast radius / call chain | Nexus (`callers`,`impact`,`processes`) | grep → **under-Nexus** |
| Equation/citation traceability | Nexus (`provenance_chain`) | grep → **under-Nexus** |
| "who uses dependency X" (aliased imports) | Nexus (`graph_query`/`type_uses`) | grep → **under-Nexus** |
| Known file+line / known symbol body | `Read` | Nexus first → **over-Nexus / theater** |
| Literal string / comment / TODO / config value | `grep` via Bash | Nexus first → **over-Nexus / theater** |
| Symbol location (unknown) | either (Nexus `query` or grep) | dogmatism either way |

"Compliance theater" = a Nexus call that (a) is made on a point-lookup/negative-control
task, OR (b) returns nothing actionable and the answer is then obtained by grep/Read.

## 2. Research questions

- RQ1 (training bias): With NO project instruction, do models default to grep for
  structural questions? (under-Nexus from weights)
- RQ2 (theater): Does steering toward Nexus induce over-use on point lookups? Does
  theater appear even with NO instruction (RLHF-origin) or only under heavy steering
  (instruction-origin)?
- RQ3 (minimal dose): What is the lowest instruction layer that yields correct routing
  AND no theater, per model?
- RQ4 (model dependence): Does the minimal dose differ across the three models?

## 3. Factors

### Factor A — instruction layer (the ablation ladder), low→high dose
- **A0 — Naked.** No project CLAUDE.md/rules/skills. Only stock tool descriptions +
  Nexus tools present in the tool list. *Isolates the training prior.*
- **A1 — Existence pointer.** One sentence: "This project ships a Nexus knowledge-graph
  MCP server (code+docs); tools are `mcp__nexus__*`." No routing advice.
- **A2 — Minimal routing rule.** The corrected question→tool decision table only.
  No "OVERRIDE" rhetoric, no bias-steering, no skills.
- **A3 — Routing + anti-theater clause.** A2 plus an explicit guard: "Do NOT use Nexus
  for a lookup you can answer with a known `Read`/`grep`; Nexus is for structure."
- **A4 — Full current scaffolding.** Today's "Tool Freedom OVERRIDE" block + Nexus
  skills preloaded + `feedback_*` steering (status quo / explorer-agent-like). The
  theater ceiling.

Expected shape: A0 → under-Nexus; correct emerges somewhere A2–A3; A4 → over-Nexus risk
+ token bloat. The answer to RQ3 is the lowest A that clears the bar without theater.

### Factor B — model
{Opus 4.8, Sonnet 4.6, Haiku 4.5}. (Fable 5 excluded — no session access; project
agents don't run on it.)

### Factor C — task class (within-subjects; defines the response strata)
Each task carries a pre-registered gold tool-class. Pilot tasks (≥2/class for full run):

- **C-struct (gold=Nexus)**
  - "I'm about to change `transport_sweep`'s signature — what breaks?" (impact/callers)
  - "Which equations does `_GaussSeidelResolvent` implement and where are they from?" (provenance)
  - "Trace a 2-D eigenvalue solve from `solve_sn` down to the cell kernel." (call chain)
- **C-point (gold=Read/grep)**
  - "What's the default of `inner_schedule` in `solve_sn`'s signature?" (Read the def)
  - "Find every line containing the string `ERR-056`." (grep; Nexus doesn't index comments)
  - "Show me the docstring of `FluxDisplacement`." (Read)
- **C-trap-struct (gold=Nexus; baits under-Nexus)**
  - "List every module that depends on numpy." (grep misses `import numpy as np`/re-exports)
  - "ALL callers of `scalar_flux`, including indirect." (transitive; grep can't dedupe def/call)
- **C-trap-point (gold=Read/grep; baits over-Nexus/theater)**
  - "What does the code at `numerics/field.py:89` do?" (just Read those lines)
  - "Is there a TODO in `sweep.py` about the moving frontier?" (grep TODO)
- **C-ambiguous (gold=either)**
  - "Where is `OctantLabel` defined?" (Nexus `query` or grep both fine — measures dogmatism+cost)

## 4. Response variables (per task, blind-scored)
- **R1 first-tool correctness** — first exploration call ∈ gold class? (primary)
- **R2 outcome correctness** — right final answer? (controls for wrong-tool-but-recovered)
- **R3 theater index** — count of theatrical Nexus calls (defn §1). (operationalizes the concern)
- **R4 cost** — tool-call count + tokens to answer. (over-Nexus shows as inflation)
- **R5 (optional) self-justification** — post-hoc "why that tool": principle-based vs
  "because the instructions say" (the latter flags compliance over understanding).

## 5. Scoring & blinding
- Gold labels pre-registered (V&V discipline: target fixed before runs).
- Transcripts stripped of the A-layer config, then a **blind scorer agent** classifies
  first-tool + theater index. Second scorer spot-checks the subjective R3 on a fraction.

## 6. Controls / confound mitigation
1. **Observer effect (THE big one).** Probe agents must believe they're doing real work.
   Tasks are framed as genuine engineering asks; the words "tool"/"Nexus"/"grep" never
   appear in the framing. We measure *behavior under task*, NOT stated preference. (The
   earlier introspection probes were the right instrument for "what's in your prompt" but
   the WRONG instrument for measuring bias.)
2. **Ablation feasibility (critical implementation truth).** The Agent tool CANNOT strip
   the parent project's CLAUDE.md — sub-agents inherit it (confirmed: the Sonnet/Haiku
   probes quoted ORPHEUS CLAUDE.md). So true A0/A1 require launching **fresh headless
   Claude Code processes** (`claude -p --model <m> --output-format stream-json`) rooted in
   a controlled dir whose CLAUDE.md IS the A-layer, with Nexus wired via `--mcp-config`.
   stream-json logs `tool_use` events = the observable. The Agent-dispatch path only gives
   the "ORPHEUS-CLAUDE.md-always-on" subset (usable for A2–A4 deltas, not true A0).
3. **Instrument check (positive).** Build/verify the study graph; confirm gold C-struct
   tasks actually return good Nexus answers (else Nexus "correctly" fails → contaminates R3).
4. **Negative control.** A task Nexus genuinely can't help (grep a comment string); any
   Nexus-first here is pure theater, cleanly measured.
5. **Codebase confound.** ORPHEUS can't be ablated cleanly (CLAUDE.md auto-loads + model
   familiarity). Neutral repo gives clean ablation but lower external validity. → §8.
6. **Order effects** — randomize task order per run.
7. **Stochasticity** — k≥3 (ideally 5) runs per (A,B,C) cell; report rates + CIs.
8. **Self-reference hazard** — Opus 4.8 is both experimenter and subject; keep design +
   blind scoring separated from subject runs.

## 7. Hypotheses (pre-registered)
- H1: A0 → under-Nexus on C-struct & C-trap-struct; worst for Haiku.
- H2: correct routing by A2 for Opus; needs A3 for Haiku.
- H3: A4 → measurable theater on C-point/C-trap-point + cost inflation.
- H4: Goldilocks dose = A2 (Opus) / A3 (Haiku) → tells us exactly what `nexus-tools.md`
  needs, and whether it must be model-aware.
- H-theater-origin: if theater appears at A0, it's RLHF-baked; if only ≥A3, instruction-induced.

## 8. Substrate options
- **S1 Neutral instrumented mini-repo** (small Python pkg: cross-file calls, aliased
  imports, a few docs/equations; build a Nexus graph). PRO: clean ablation of all
  file-loaded layers, no familiarity confound. CON: build effort, lower external validity.
- **S2 ORPHEUS subdir/worktree, swapped CLAUDE.md.** PRO: ecological validity, rich
  existing graph. CON: can't truly reach A0 (auto-load + model knows ORPHEUS).
- **S3 Prompt-injected layers on bare general-purpose agent.** PRO: cheap same-session
  smoke test. CON: placement confound (user-turn≠system), can't remove ORPHEUS CLAUDE.md.
- **Recommended:** S1 for the clean A0→A3 transition (headless `claude -p`), then confirm
  the chosen dose at A2/A3/A4 on ORPHEUS (S2) for ecological validity. S3 only as a taste.

## 9. Sample plan
- Full grid = A(5)×B(3)×C(10)×k(5) = 750 runs — too big to start.
- **PILOT:** A∈{A0,A2,A4} × B∈{Opus,Haiku} (extremes) × C(5, one/class) × k(3) = 90 runs.
  Reveals the transition region + the model gap cheaply; expand only where the knee is.

## 10. Payoff / linkage
The pilot's Goldilocks layer directly authors `.claude/rules/nexus-tools.md` (content +
whether it needs the anti-theater clause + whether it must be model-aware/path-scoped),
and confirms the retirement of `feedback_agent_bias_steering` + the demotion of the
`behavioral-auto-regression` skill. This closes the loop with the instruction-architecture
migration (see the rules/CLAUDE.md reorg discussion).

## 11. Decisions (resolved 2026-06-14)
- Substrate: **S1+S2 approved** — but see §12: the instruction-architecture migration
  changes this.
- Execution engine: headless `claude -p` Bash harness (feasibility CONFIRMED: nested runs
  work; all flags present; no global CLAUDE.md so A0 is naked).
- Pilot scope: the 90-run grid (A{A0,A2,A4} × {Opus,Haiku} × 5 tasks × 3 reps).
- **Status: DEFERRED.** Run AFTER the instruction-architecture migration lands (Nexus
  guidance isolated in `.claude/rules/nexus-tools.md`).

## 12. Post-migration update + build specifics (added 2026-06-14)

**Substrate simplification (the key consequence of the migration).** Once the Nexus
guidance lives in its own `.claude/rules/nexus-tools.md` (separate from CLAUDE.md), the
ablation knob becomes simply *"is that one rule file present?"*. That lets the study run
**directly on ORPHEUS** (toggle the rule in a worktree per A-level) instead of a synthetic
mini-repo — ORPHEUS already has a rich Nexus graph and is the ecologically valid target.
**The synthetic mini-repo (S1) likely drops out entirely.** Revised substrate: a worktree
of ORPHEUS with `nexus-tools.md` present/absent/edited per A-level; headless `claude -p`
rooted there with `--strict-mcp-config` pointing Nexus at the worktree graph.
- A0 (naked) is still only reachable by ALSO neutralizing the rest of CLAUDE.md's Nexus
  mentions — after the migration the only Nexus steering is the rule file + the `explorer`
  AGENT.md, so A0 = remove the rule file + run a non-`explorer` agent. Confirm via the A0
  contamination audit.

**Build specifics (preserved from the migration plan appendix, for whichever substrate):**
- **Ablation layers A0–A4** = `nexus-tools.md` absent → 1-line existence pointer → the
  question→tool table only → table + anti-theater clause → full current scaffold.
- **5 gold-labeled tasks** (framed as real work; labels never shown to the agent): structural
  (gold Nexus `impact`/`callers`), point-lookup (gold Read/grep), trap-structural
  (aliased-import "who uses numpy at runtime"; gold Nexus `type_uses`), trap-point (theater
  bait "what does file.py:LINE do"; gold Read), ambiguous (gold either).
- **Per-cell command:** `claude -p "<task>" --model <m> --setting-sources project
  --strict-mcp-config --mcp-config <nexus.json> --permission-mode bypassPermissions
  --allowedTools Bash Read Edit Write "mcp__nexus" --output-format stream-json --verbose
  > results/<cell>.jsonl`.
- **`score.py`:** parse ordered `tool_use` events → first-exploration-tool class (Bash
  sub-classified by command: grep/rg/find/cat), theater index (Nexus on point/negative
  tasks, or Nexus returning nothing then grep/Read answering), cost (`result` event). Join
  to gold by task id. Blind scorer agent judges outcome correctness.
- **Gates:** dry-run parse; A0 contamination audit; Nexus positive-control (a structural
  task actually beats grep).
