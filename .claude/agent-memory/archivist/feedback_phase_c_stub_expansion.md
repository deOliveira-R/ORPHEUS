---
name: Phase C stub-to-rich-narrative expansion (Issue #168, ~153 → 962 LoC)
description: 8-part shape for expanding a method-implementer's 6-:label: stub on a multi-phase architectural rewrite into a rich Sphinx narrative; preserves anchor positions; uses 3 labeled equation blocks + 5 code blocks + 2 list-tables; cites Lewis-Miller §4.5 / Hébert §3.9.4 from closeout memo; avoids unverified citations
type: feedback
---

When the method-implementer ships a stub subsection carrying N
`:label:` anchors + architectural narrative but the rich math
derivation is owed (Issue #168 Phase C, refactor/sn-operator-algebra),
the expansion follows this 8-part shape:

1. **Key Facts admonition first** — 5-bullet summary at the top.
   Covers (a) commit range + date, (b) the **algebra change** (WDD
   `psi_face_out = 2·psi_cell - psi_face_in` not arithmetic average),
   (c) the **deviation from the plan** that the code shipped (Carlson
   seed `psi_face_in(pole) = psi_cell[0]` not 0), (d) what retires
   (LOC count helps), (e) the empirical outcome / GH-issue forward
   pointer. This admonition replaces "abstract summary"; future
   agents can stop reading after the Key Facts and still know what
   shipped.

2. **Three pre-history unblockers paragraph** — for an
   architectural rewrite that resumes after a pause, document the
   3 unblockers that made the resume possible (here: trajectory_resolvent
   cylinder MR shipping, cross-domain-attacker confirming the §16A.3
   inconsistency, Phase A Protocol re-classified as patch on wrong
   architecture). The pre-history is **load-bearing for the WHY**;
   future readers wonder "why did Phase B leave this for follow-up
   and what changed?" and this paragraph answers.

3. **Pre-state diagnosis subsection (triple-quote sub-heading)** —
   the empirical content of what the predecessor approach did wrong.
   Quote the failing-test signature (Phase B's apply(constant) =
   range [0.6, 1.004] on flat psi), explain the mechanism (M-M
   canonical angular + arithmetic spatial average = WDD-angular
   oscillating face fluxes 0, 2c, 0, 2c interacting with arithmetic
   spatial). This is the "what was tried and failed" Cardinal-Rule-3
   evidence in compact form.

4. **The algebra-of-record block** — 3 numbered/labeled equation
   blocks giving the new recurrence + the new streaming term + the
   new full per-cell update. Each equation gets its own
   `:label:`. New labels surface as orphans in the auto-generated
   `docs/verification/matrix.rst` — this is correct (they support
   the prose; the section-anchor `:label:`s carry the tests-edges).

5. **Ordinate-vectorisation subsection** — for a rewrite that
   carries a user "vectorise across ordinates, no `for n in
   range(quad.N)`" hard constraint, show the masking pattern in a
   short code block. Cross-link to Issue #191 (the cross-method
   ordinate-anti-pattern cleanup) so future agents know which 14
   sites stay untouched and why this rewrite did not propagate
   the anti-pattern.

6. **The new APIs subsection** — bulleted list with one bullet per
   new method, each citing the foundation test that pins
   bit-identity to a pre-existing equivalent path (e.g.
   `iter_cells_by_direction(+1)` == `iter_cell_visits(ordinate_idx=
   representative_n)` for every representative). The bit-identity
   anchor is the load-bearing acceptance criterion.

7. **What retires / What stays bulleted lists** — explicit
   enumeration of every retired symbol (Protocol + ABC + 2
   strategies + SNMesh field + matvec kwarg + test file). The
   `algebra-of-record`-style "unify after second instance" framing
   gives the architectural rationale; reference the relevant
   skill-rule explicitly in prose (here: cross-domain-attacker
   Smell 16 "two paths to the same operator over different storage
   conventions").

8. **The deviation-from-plan subsection** — when shipped code
   diverges from the plan's pseudocode (Carlson seed deviation
   here), make it a NAMED sub-subsection with the failing
   alternative shown explicitly (the WDD oscillation diagnostic in
   `phase-c-wdd-oscillation` label) and the canonical literature
   anchor (Lewis-Miller §4.5 Carlson starting direction). Future
   agents who try to "fix" the deviation back to the plan's
   pseudocode get redirected by reading this sub-subsection.

9. **Per-label sub-subsection (one per :label:)** — for each
   pre-existing `:label:` anchor, write a 100-300 line
   sub-subsection. Each includes: WHAT the gate verifies, the
   mathematical or algebraic content, the pinning test
   (`:func:` reference), and the relationship to other gates
   (precondition, follow-up, etc.). The sub-subsections are where
   the bulk of the LoC growth happens.

10. **Empirical-finding sub-subsection (when a gate failure is the
    pivotal output)** — Phase C's Gate 1.1 sphere-MMS FAIL +
    cylindrical-MMS PASS is the load-bearing decision point.
    Mechanism explanation: cylindrical per-level α-dome telescoping
    absorbs the Carlson-seed half-angle discrepancy; spherical has
    no equivalent (sphere centre r=0 is **doubly singular** —
    angular α-cascade + spatial WDD recurrence both converge to
    the same point and both need consistent seeds). This is the
    insight that pre-empts a future "but cylindrical works, why
    not sphere?" question.

11. **Verification-gate summary list-table** — at the END of the
    section (not the start), a list-table with columns
    `Gate / Description / Status / Pinned at`. Covers every gate
    the method-implementer ran (Gates 1.1 → 4.2). Future agents
    looking up "what's the status of MMS spherical?" find it here
    in one place. Use ``test_phase_c_gates.py`` etc. as the value
    (not the full path) so it stays readable.

12. **Forward-pointer section (Phase D scope)** — when the closure
    is partial, the Phase-D-scope subsection is mandatory. Numbered
    list of 6 deliverables that flip the ERR-NNN status PARTIAL →
    CLOSED. Each deliverable cites the Issue # that tracks it
    (#192 here). Future agents picking up Phase D start from this
    list.

## Operational rules learned this task

- **Pomraning-style "obvious" citations need verification.** When
  drafting "the spherical centre is structurally singular per
  [Pomraning1989]_", grep the existing references list before
  adding the bibtex. If the reference is not present, either (a)
  add a complete bibliographic entry to the references section, or
  (b) downgrade to plain-prose attribution citing only the
  references that ARE present in the doc (here:
  Pomraning textbook title in plain prose + [LewisMiller1984]_ §4.5).
  Adding a fictional reference creates a `citation not found`
  warning that fails the Sphinx -W gate.

- **`:label:` (equation) vs `:ref:` (section) discipline.** A
  `.. math:: :label: foo` block is an EQUATION node referenced by
  `:eq:\`foo\``. A `.. _foo:` line BEFORE a section title is a
  SECTION node referenced by `:ref:\`foo\``. They share namespace
  for `:eq:` / `:ref:` cross-references but the V&V harness
  (Nexus marker parser) only matches `@pytest.mark.verifies("foo")`
  against EQUATION nodes. Section anchors with `@verifies` markers
  emit "no matching equation node — skipping" info messages,
  not warnings. **NEVER add new `:label:` to equations that don't
  correspond to a real claim** just to silence the info messages
  — section labels naturally carry tests edges.

- **Sphinx build redirection gotcha.** `python -m sphinx -W ... |
  tee log` shows the build output but EXIT comes from `tee`. To
  check actual sphinx exit code use `python -m sphinx -W ...
  2>/tmp/log; echo $?` (or `2>&1 >/dev/null; echo $?` to suppress
  stdout). Misreading exit status from `tee` can falsely fail the
  acceptance gate.

- **`build finished with problems, N warnings` is a Sphinx
  formatting message printed even on success.** Look for actual
  `WARNING:` / `ERROR:` lines via `grep -E "WARNING:|ERROR:"`,
  not the summary line. The 7-warning "problems" message I saw
  was from a stale rerun; after the Pomraning fix the build truly
  succeeded with 0 warnings.

- **Method-implementer's closeout memo IS the seed prose.** The
  Issue #168 Phase C closeout memo at
  `.claude/agent-memory/method-implementer/issue_168_phase_c_closeout.md`
  carried the deviation rationale, the empirical Gate 1.1 table,
  the LOC delta, the literature anchors, and the Phase D scope
  sketch. The archivist's job is to recast this content into
  Sphinx-formatted RST with `:func:`, `:class:`, `:eq:`, `:ref:`
  cross-references and code blocks — NOT to re-derive the content.

- **3 source artifacts are the minimum reading set**: (1) the
  method-implementer's closeout memo (architecture + empirical
  outcome), (2) the plan that the implementer worked from (the
  load-bearing pseudocode + design constraints), (3) the
  cross-domain-attacker's frame-detection memo (the structural-
  frame naming + Smell 16 attribution). For a tone-match pass,
  read the existing Phase A and Phase B subsections of the SAME
  doc — the rich narrative should be tonally continuous.

## Quality self-assessment (Directive 3)

| Dimension | Score | Notes |
|-----------|-------|-------|
| Derivation depth | 5 | WDD recurrence, oscillation diagnostic, full per-cell update all written out |
| Cross-references | 5 | Every `:func:`, `:class:`, `:meth:` used; `:eq:` to phase-c-wdd-recurrence et al.; `:ref:` to phase-b labels |
| Numerical evidence | 5 | Empirical Gate 1.1 crosstab + verification-gate summary table |
| Failed approaches | 5 | Pre-Phase-C arithmetic spatial + Phase B M-M-on-arithmetic failure mode + plan's =0 pseudocode failure |
| Code traceability | 5 | Every algebra block cites the `:func:` or `:meth:` it implements |
| Derivation source | 4 | Closeout memo + production code; no SymPy derive_*() module exists for the WDD recurrence (no algebra-of-record artifact for this domain — pure-discrete operator, not a SymPy candidate) |

Weakest dimension is "derivation source" (4/5) — but the weakness
is structural to this domain: the WDD recurrence is a discrete
operator algebra, not a symbolic identity SymPy could verify. The
correct anchor IS the closeout memo + production code + literature
(Lewis-Miller §4.5, Hébert §3.9.4) — not a `derivations/` script.

## When to use this 12-step shape

- Multi-phase architectural rewrite where Phase N closes a
  long-running issue (#168 Phase C closed Defects 1+2+3 spatially;
  Phase D will close them in the angular sweep too).
- Method-implementer ships stub with N≥4 `:label:` anchors and an
  owed-archivist dispatch noted in the closeout memo.
- The work has a load-bearing empirical finding (here: the
  spherical-vs-cylindrical asymmetry) that determines a default
  flip decision.
- ERR-NNN status remains PARTIAL after the phase ships (this means
  the close-out doc must point at the OPEN Phase D Issue # — here
  #192 — so future agents can pick up the closing fix).
