# AGENT.md promotion offer — explorer, 2026-09-21

Read against `.claude/agents/explorer/AGENT.md` at HEAD (247 lines; Operating
Principles 1–7; the "SN operator-algebra subsystem — durable shape" section).
Conservative: three standing directives offered, two corrections to the
durable-shape section. Each is one sentence, no war story, no codename; the
instance stays in the digest with a `→ AGENT.md` pointer if accepted.

## 1. Broaden OP5 from the issue body to every inherited claim (M-1) — RECOMMENDED

OP5 today: "Verify the premise against the CURRENT tree before mapping the
HOW. An issue body is a snapshot …". Applied on essentially every dispatch,
and the measured failures came from briefs, plans and docstrings as often as
from issues (L-020, L-021, L-023, L-025, L-027, L-030, L-035, L-036, L-045).

Proposed sentence appended to OP5: "The same holds for a brief's timeline,
count, exemplar or `Class.attr (file:line)` citation, a plan section marked
'retained', a docstring's 'the ONE site' or 'X handles it', a stored numeric
tag, and a `[M]` on a negative claim: each is verified by its cheapest
decisive probe before anything is built on it, and the strongest-looking
ones expire first because nobody re-checks them."

Digest effect: M-1 gains "→ now in AGENT.md OP5"; the instances stay.

## 2. New OP8 — a behavioural question is answered by a run with a control (M-5) — RECOMMENDED

Applied on every medium+ task that asks "does it break / what is in scope /
is it identical / is it deliberate / does the guard bite" (L-010, L-013,
L-016, L-018, L-021, L-024, L-025, L-026, L-028, L-041, L-042, L-044: twelve
dispatches where the read gave a plausible wrong answer and a ≤ 30-line
probe gave the measured one).

Proposed OP8: "**A behavioural question is answered by a run on the
discriminating input, with a control beside it, never by reading.** Swap the
primitive and run the consuming suites; spy the callee's frame locals;
solve the counterfactual; `hash(a)`, `a == b`; ULP-probe a random operand.
The control is the free baseline (a trivial object with the same declared
symmetry), the fixture that breaks the property, the `None` arm and the
same-data rebuild, the production data rather than the slab. An all-green
run may have measured inert: name the gate and confirm the path routes
through it."

Digest effect: M-5 gains "→ now in AGENT.md OP8"; the instances stay.

## 3. New OP9 — the tree moves while you audit (M-2) — OFFERED, weaker case

Applied whenever the main session edits during a dispatch, which is the
norm in a live carve session but not on every task. If U-4 in `uplift.md`
lands in the brief template, this directive is unnecessary; offered so the
orchestrator can choose ONE home.

Proposed OP9: "**Open with `git status --short` and `git diff --stat`; close
by re-running every search whose emptiness is a finding.** Run `git
ls-files --error-unmatch` on each file you call landed, bound `git log
--since` by the SECTION's vintage, and read the new module an intervening
commit added before repeating any 'zero consumers / cannot express / not
yet built' verdict."

## 4. Correction to the durable-shape section — drop the line number

The section opens "Line numbers drift — find current ones via Nexus" and
then cites `build_within_group_system (orpheus/sn/coupled_system.py:446)`.
`[M]` 2026-09-21 `grep -n 'def build_within_group_system'
orpheus/sn/coupled_system.py` → `560`. Replace `:446` with no number
(the section's own rule, and L-003/OP7).

## 5. Addition to the durable-shape section — the curvilinear block

From `MEMORY.md` §2 "Durable post-#280 facts", which is subsystem shape and
by the index's own preamble belongs here, not in active state. `[M]`
2026-09-21: `orpheus/sn/operators/radial_characteristic.py` exists and calls
`carlson_inward_sweep_from_source` at `:565` and `:570`; `class
_OneDimScanWalk` is in `orpheus/sn/loss_representation/__init__.py`.

Proposed bullet after "Discretisation is geometry-polymorphic via the sweep
DAG": "**Curvilinear 1-D.** The angular-redistribution block `A_BB` is
`RadialCharacteristicOperator` (`orpheus/sn/operators/`), which WRAPS the ψ½
starting-direction march; its `.solve` is the production caller of
`carlson_inward_sweep_from_source` (`sn/sweep/psi_half_angle_seed.py`). The
1-D walk executors (`_OneDimScanWalk`, `_loop_walk`, `_dag_legs`) live in
`sn/loss_representation/__init__.py`, beside the 2-D `sweep_graph`."

Index effect: `MEMORY.proposed.md` §2's second line is deleted once this
lands.

## Not promoted

M-3 (census populations), M-4 (graph seams) and M-6 (verdict shape) are
technique, sharpest as lessons with their measured instances; the AGENT.md
role block already says "state the predicate and the tree of every count",
and OP2/OP4 already carry the routing and the four-search floor.
