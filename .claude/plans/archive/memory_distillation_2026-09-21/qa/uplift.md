# Uplift candidates — corrections that generalise past this agent

Five. Each names the rule or skill it would join, the clause in that file's own
form, and the founding case. The first four were already carried as "standing
debt" in `qa/MEMORY.md` and are re-verified here by READING the current files
(2026-09-21): none of the five is present. The four that HAD landed since the
debt table was written (`#17(a)` per-consumer-kind, `#17(i)` prescribed repair,
`#17(g)` control-per-stage, `#34` α-normalised AST, `#35` overloaded unit,
`#36` deselected catcher) are struck from the debt and now appear as §I rows in
the digest.

---

## U1 — a clause naming TWO mechanisms under ONE check

**Home:** `vv-principles` § Anti-patterns (a new numbered item), or
`instrument-doctrine` X3, which already owns "prose is not enforcement".

> **NEVER** accept a clause, marker or docstring that names TWO mechanisms under
> ONE mechanical check — the check reaches one and the prose credits both.
> **check:** for every mechanism the text names, point at the line of the check
> that would find it. **tell:** an imperative that says "X plus Y" whose check
> mentions only X.

**Founding case:** a §6b census spelling named "a shape minted independently by
consumers **plus** `isinstance` doors on the producer's type" under the single
check "grep the shape's CONSTRUCTOR", which structurally cannot find the doors;
the source had stated them as *two* spellings (L-081; digest E10). Verified
absent: `vv#25`'s "two mechanisms = two checks" is about a bundled change's NULL
result, a different subject.

## U2 — a retirement or absorption NOTE is a CARRIER claim

**Home:** `retirement-audit` § E (retiring the retirement's own residue), beside
item 19, or `vv-principles` § Anti-patterns.

> **NEVER** accept a retirement note ("absorbed into X", "superseded by Y") on
> the strength of the TITLE match — a lesson, rule or clause is a BUNDLE, and the
> note names the one mechanism its title shares with the carrier.
> **check:** count the SOURCE's own rules, then grep the carrier for the most
> distinctive token of EACH; read the source's closing paragraph, which is where
> a body states what still holds. **tell:** a note written from the title; a
> "mechanism superseded" verdict on a body whose last line says why it stays.

**Founding case:** 20 retired lessons audited — 12 carried FULL, 8 PARTIAL, ~14
mechanisms unstated, one (a commit-message backtick-substitution hazard) with 0
hits tree-wide after its "absorbed" note (L-083; digest A23).

## U3 — a recall counter DOWNSTREAM of a filter, and an inferred relation under a declared name

**Home:** `instrument-doctrine` X1 (beside the input-count clause landed
2026-09-21) and X2 (the census protocol's decoder step).

> A counter placed downstream of a FILTER cannot count what the filter dropped:
> a total of zero is compatible with "nothing matched" and with "everything was
> discarded before the counter".
> **check:** require a per-REASON drop breakdown between input and output, and
> assert the input count separately. **tell:** `found: 0 / unresolved: 0`,
> exit 0, on an artifact you know is non-empty.
> Companion: a relation whose rows are INFERRED (a name-token guess) published
> under a DECLARED relation's name — read the predicate that sets the status
> before quoting the number.

**Founding case:** `nexus runtime-ingest` printed `nodes: 0 / unresolved: 0` on
a real coverage report whose 339 file keys were all dropped by an
absolute-vs-relative path filter; normalised, the same artifact joined 2892
nodes (L-070; digest A11). The inferred half: all 16624 `implements` edges are
`source="inferred"`, 81 % on one shared token, while the status reads
"verified" (L-070; digest E7). The X1 clause that landed today covers an EMPTY
INPUT LIST, not a filter that silently empties a non-empty one.

## U4 — a before/after `[M]` pair states ONE instrument for both halves

**Home:** `plan-authoring` § 4 (the markers section), as the temporal twin of
RATIO-NEEDS-ITS-POPULATIONS.

> **BEFORE/AFTER-ONE-INSTRUMENT** A before/after pair is a ratio in time: both
> halves are read with the SAME instrument, named beside them.
> **check:** name the command that produced each half; if they differ, there is
> no delta. **tell:** a "≈29.5K → ≈14K" pair whose second half is reproducible
> from a recorded tool run and whose first half matches no recorded instrument.

**Founding case:** `rules/workflows.md:70`'s `[M]` "≈29.5K" matched no recorded
instrument while its "≈71K" half was exactly keep−omit from a different one
(L-086; digest E16).

## U5 — a `.claude/` sweep belongs in a retirement's blast radius

**Home:** `retirement-audit` § B (surfaces a symbol grep cannot reach), as a new
item beside B.4.

> A type or concept retirement's blast radius includes `.claude/agents/*/AGENT.md`,
> `.claude/skills/*/` and `.claude/agent-memory/*/`. An AGENT.md outranks a
> production docstring: it loads FRESH on every dispatch, so a stale brief is
> re-injected as current fact and the sub-agent's output is indistinguishable
> from a correct one.
> **check:** grep the retiring name across `.claude/` with the code grep, and
> sort the hits by TENSE (item 19). **tell:** a retirement whose doc sweep
> covered `docs/` and `orpheus/` only.

**Founding case:** CS3-R — 3 of 12 surviving references to a retired type were
agent briefs (one teaching a retired 4-role grid, one carrying as an imperative
what its own source had already ⛔-corrected). Agent memory was the biggest and
least-swept slice: 182 lines / ~20 files, against 75 lines for
skills + agents + rules combined (L-071; digest D15).
