# Archivist topic-file blast-radius audit — proposal

**Dispatch:** the harness campaign's memory-distillation close-out, owed item (1). Read-only on
tracked files; this file is the only deliverable. **28 candidates, each read in full.**

**Tally: 0 KEEP · 15 RETIRE · 13 SALVAGE-then-RETIRE.**

## The headline finding

The 2026-09-21 distillation judged all 28 "archaeology" **by name, without reading them**. Read,
the by-name judgement is right about the *file* in 28 of 28 — no candidate is a live pointer, a
standing convention or a load-bearing referent — but **wrong about the contents in 13 of 28**,
which carry something nowhere else in the corpus. Two of those 13 (`425_sn_chapter.md`,
`425_outside_chapter.md`) are the **only** record of an entire pass: `lessons_archive.md` runs
L-001…L-113 and has **no section for #425** (it jumps L-099 → L-103 over 2026-09-07). Retiring
them unread would have lost it.

Three claims carried by these files are **false today**, and carrying them forward is a net
negative — that is an independent argument for retirement, not merely an absence of value. They
are listed under "Refuted claims" below.

---

## Method, and what was verified

**Every status line was verified against git / GitHub, never against a memory's claim.**
`[M]` 2026-09-21: `9ac5269e`, `13bcf5f`, `6d16df6`, `3356cec`, `87a3fed` are all ancestors of
`HEAD`; branch `docs/425-within-group-algebra` is gone; `gh issue view` reads **CLOSED** for
#425, #196, #195, #247, #251, #257, #168 and **OPEN** for #240 and #263 (see the note on the two
open ones below); `git status --porcelain -- docs/` is empty, so every page cited below is
committed.

**The referrer census was re-run independently, with a positive control.** `[M]` 2026-09-21, a
python scan (`re`-free substring match over every text file of the repo minus
`.git`/`_build`/`.venv`/`scratch`, plus the main agent's memory dir): **2 430 files read**,
positive control `feedback_canonical_convention_page` → **2** referrers (non-zero, so the filter
works). The result **exactly reproduces `referrers.md`** — the same 11 non-table lines, no more
and no fewer — and each of those 11 was additionally confirmed by re-running its own
`grep -n "<stem>" <file>`. The two censuses' denominators differ (2 430 vs the brief's 2 242)
because the predicates differ — mine reads the main-agent memory directory and excludes more
binary suffixes — not because they disagree; the **finding** is identical under both.

**`MEMORY.md` needs no edits at all.** `[M]` none of the 28 stems appears in
`.claude/agent-memory/archivist/MEMORY.md`. The candidates were never in the hot index, which is
consistent with the distillation reaching them by directory listing rather than by index.

**Three measurements were taken for this audit**, because each decides a salvage and each would
otherwise have entered the digest's §2a "what warns / what is silent" list as a relayed claim —
the exact failure L-113 exists to prevent. All are throwaway one-page Sphinx projects in this
venv, each with a two-sided control.

| # | claim | verdict |
|---|---|---|
| M1 | a citation defined and **not** referenced | **WARNS** at default severity: `WARNING: Citation [Orphaned1999] is not referenced. [ref.citation]`; the referenced citation on the same page was silent |
| M2 | a bare `Γ_-` / `S_-` in prose | **ERRORS** at default severity: `ERROR: Unknown target name: "γ"` / `"s"`; the escaped `Γ\_-` and the no-trailing-underscore `V_bulk` control were both silent. `-W` **exits 1** on it (clean-control page exits 0) |
| M3 | an **over-long** section underline | **silent**; only *too-short* warns (`Title underline too short`), both observed in one build |

---

## Verdict table

`table-only` abbreviates: the file's sole inbound line is
`.claude/plans/archive/memory_distillation_2026-09-21/archivist/table.md`, the distillation's own
candidate list — a frozen, already-past-tense historical artefact that stays true whatever happens
to the file, so it needs **no edit**. Referrer edits are numbered **R1–R5** and specified in full
below; there are only five, because every other inbound line is itself inside a retiring
candidate and dies with it.

| file | verdict | reason (one line) | referrer edits |
|---|---|---|---|
| `feedback_plan_triage.md` | RETIRE | Its whole content — the 5-action DELETE/UPDATE/POST/CONSOLIDATE/SUPERSEDED rubric, archive-vs-`rm`, the master-research-log exception, the status-update-vs-plan-summary anti-pattern and both worked clusters — is the `plan-cluster-triage` skill (§"The 5-action rubric", §"Defaults", §"Status-update vs plan-summary anti-pattern", §"Worked example"). | table-only |
| `feedback_misc_cluster_triage.md` | RETIRE | Its four sub-patterns, the `gh api -X PATCH` recipe, the "issue OPEN despite Closes-#NNN" pattern and the citation-hygiene header are `plan-cluster-triage` §"The 4 sub-patterns for heterogeneous batches", §"The PATCH existing comment sub-pattern", §"The issue OPEN despite Closes-#NNN comment pattern" and §"Citation hygiene for relocated content". | table-only |
| `feedback_plan_triage_quadrature_cluster.md` | RETIRE | The quadrature cluster IS the skill's §"Worked example — post-#138 quadrature cluster (2026-04-30)", including the 4-of-4 + 2-of-2 DELETE count, the precision-floor exception and §"Asymmetry rule"; §"Cross-check before deletion" carries its closing grep. | table-only |
| `feedback_phase2a_large_relocation.md` | RETIRE | The ≥900-LoC anchor-preserving stub, its five steps and its counter-pattern are `doc-issue-relocation` §"The ≥ 900-LoC stub pattern" Steps 1–5 + Pitfall 12, verbatim down to the stub skeleton. | table-only |
| `feedback_phase2b_label_dense_partition.md` | RETIRE | Its 4-class anchor triage is `doc-issue-relocation` Step 1 + Step 5 (class 3, "matrix.rst-only ⟹ droppable", is that skill's own closing paragraph); its stale-text-flip is digest §4; its commit cadence is superseded by digest §9's build-TWICE sequencing; and its counter-pattern carries **refuted claim F1**. | table-only |
| `feedback_phase2c_staleness_sweeps.md` | SALVAGE → RETIRE | Salvage **S1**; the rest is `doc-issue-relocation` Pitfall 8, digest §4 (capability flip / forward-framing) and §9 (`grep -v _build`), and its "`-W` after every group" cadence is superseded by digest §9. | table-only (its one other inbound is `feedback_phase_d_carlson_seed_narrative.md:153`, itself retiring) |
| `feedback_phase2d_table_preservation_constraint.md` | SALVAGE → RETIRE | Salvage **S2**; table preservation itself is `doc-issue-relocation` Pitfall 8 and keep-the-anchor-through-a-retitle is digest §4. | table-only |
| `425_sn_chapter.md` | SALVAGE → RETIRE | **#425 has no `lessons_archive.md` section** — these two files are the pass's only record. Salvage **S3–S5** into a new archive section **L-114**. `[M]` `9ac5269e` is an ancestor, the branch is gone, #425 is CLOSED. | table-only |
| `425_outside_chapter.md` | SALVAGE → RETIRE | The richer half of the same pass: salvage **S6–S9** into L-114. Its ⛔ item 2 is **refuted claim F2** and the refutation is itself the salvage (**S7**). | table-only |
| `feedback_bc_trace_law_wave_12.md` | SALVAGE → RETIRE | Salvage **S10** (M2); the rest is Cardinal Rule 3, digest §3 (label namespaces), §6 (zero-consumer mint) and §9 (count-unchanged, since sharpened to SET-unchanged); its "build incrementally" advice is superseded by digest §9; its Pitfall 1 is **refuted claim F3**. | **R1**; two other inbound lines are inside retiring candidates |
| `feedback_wave_o_operator_algebra_docs.md` | SALVAGE → RETIRE | Salvage **S12**, **S13**; its breadcrumb-resolution half is digest §4 in sharper, census-carrying form, and its per-step kernels are the content of `operator_algebra.rst` — `[M]` the family anchor `.. _bc-extraction:` is live at `docs/theory/foundations/boundary_conditions.rst:2780`, so the page is the record (Cardinal Rule 3). | table-only |
| `feedback_post_wave_cleanup_docs.md` | RETIRE | The event class is digest §6's REFUSAL-BECOMING-A-CAPABILITY five moves (L-076) and §4's capability-flip / seam-discharge clauses, both later and sharper; its one original note (the arc's "novel extensions falsified" step is *absent* here, say so) is AGENT.md's arc preamble. | **R2**; its other inbound line is inside a retiring candidate |
| `feedback_err058_success_closeout_supersedes_phase_chain.md` | SALVAGE → RETIRE | Salvage **S14**; the arc is digest §6's fix-by-retiring SUCCESS-RESOLUTION chapter (L-013), and ⚠ its `audit --strict` baselining recipe prescribes **`git stash`**, which `process-discipline` forbids and digest §9 replaces with `git archive HEAD` into a temp tree. | **R3** |
| `feedback_issue_196_eigenvalue_verification_closeout.md` | RETIRE | Its one durable ruling is already durable in the corpus: `[M]` `docs/theory/methods/sn/curvilinear_numerics.rst:2995` carries `.. _sn-issue-196-bit-identical-vs-floor:`, cited from `:2139` and `:62` and echoed at `error_catalog.rst:1796`. Everything else is digest §4 and §7. Raised as skill-uplift **U1**. | **R4** |
| `feedback_issue_257_consolidated_taxonomy_pass.md` | RETIRE | Memo-names-a-symbol-that-does-not-exist is digest §1a; vv-status on structural labels is §3; "`-W` catches a dangling intra-doc `:ref:`, so place the anchor in the same edit" is §2a (now known to hold cross-doc too); the `operator_algebra.rst` marker ladder is in AGENT.md verbatim. | table-only |
| `feedback_d5b_d6_campaign_closeout.md` | SALVAGE → RETIRE | Salvage **S16** (M3); the ERR-061 content is `error_catalog.rst` (ERR-061, 3 catchers) plus the page, the eq-label/orphan policy is digest §3, and necessary-not-sufficient-matvec is `vv-principles` #5/#37. **#240 is OPEN** — see the note below. | table-only |
| `feedback_issue_247_legA_mode10_closeout.md` | SALVAGE → RETIRE | Salvage **S17**, the diagnostic half of the marker-ladder rule AGENT.md lacks; the Mode-10 content is `vv-principles` #10. | table-only |
| `feedback_issue_251_legB_boundary_trace_closeout.md` | RETIRE | **Its own owed action is discharged.** LESSON 1 says *"propose it to vv-principles"*; `[M]` that clause is live at `.claude/skills/vv-principles/SKILL.md:230`. LESSON 2 is digest §6's L-005 read-order; LESSON 3 is S17's positive control. | table-only |
| `feedback_issue_257_s9_coherent_promise_criterion_seam.md` | SALVAGE → RETIRE | Salvage **S19**; its homing lesson is digest §6 three times over (L-079 / L-064 / L-057); `[M]` its content is live at `docs/theory/foundations/operator_algebra.rst:5529` and `docs/theory/verification/sn.rst:3169`; its rider *"cross-doc dangling renders plain-text, `-W` blind"* is refuted by digest §2a. **#263 is OPEN** — see the note below. | table-only |
| `feedback_phase_c_stub_expansion.md` | SALVAGE → RETIRE | Salvage **S18**; the 12-part stub shape is digest §6's L-005 read-order, citation-before-adding is §9, and "do NOT mint an eq-label to silence a `no matching equation node` INFO line" is AGENT.md verbatim. | table-only |
| `feedback_phase_d_carlson_seed_narrative.md` | SALVAGE → RETIRE | Salvage **S20**; ⚠ its "Equation-label discipline" section states a **false mechanism** (refuted claim F4). | **R5** (cross-owner); its other inbound line is inside a retiring candidate |
| `feedback_s10a_emission_spectrum_value_object.md` | RETIRE | `[M]` the work SHIPPED and is committed: `docs/theory/foundations/cross_section_data.rst` carries `.. _emission-spectrum-simplex-law:` (`:1959`), `emission-spectrum-fission-source` with its `vv-status` (`:2148`/`:2152`) and `emission-spectrum-chi-mix` (`:2322`), and **0** `Archivist expansion needed` remain — so the memo's "NOT committed" line is stale and the page is the record. Its vv-status / foundation-not-L-level rulings are digest §3 and §7. | table-only |
| `issue_196_pr_index_6_closeout.md` | RETIRE | A 708-line session report (diffstat, test paste-backs, a 12-row criteria table, a session trail); its one durable contribution, the §8.1 keep-vs-flip rubric, is already distilled as a 5-row table in the **surviving** `feedback_canonical_convention_page.md:101-109`, together with the page anatomy. Its "Next step pointer" (PR-INDEX-7) was declared moot by `task57`. | table-only |
| `issue_196_pr_cleanup_docs_closeout.md` | RETIRE | The same shape one PR later — a gate audit, an anti-recommendations acknowledgement and a session trail; the phantom-label resolution discipline is digest §3, and the campaign's durable form is `feedback_canonical_convention_page.md`. | table-only |
| `peierls_greens_phase123_rich_narrative.md` | SALVAGE → RETIRE | `[M]` the narrative landed (`docs/theory/references/trajectory_resolvent.rst`, 19 hits on the family anchors) and ERR-035 is catalogued at `error_catalog.rst:2433`. Of its three proposed `vv-principles` uplifts, **two are already covered** (uniform-source blindness by Mode 7's activates/nulls declaration; convergence-rate-as-fingerprint by `numerical-bug-signatures` Signatures 6 and 7, in sharper discriminator form) and one is not — raised as **U2**. | table-only |
| `task55_transport_sweep_docs_pass.md` | RETIRE | A 56-site disposition table; its three dispositions are digest §4's finer five-register tense sort, and its MINT-vs-REPOINT adjudication for a phantom `verifies` target is digest §3's L-003 clause. Its one owed follow-up (`psi_bc` needs its own pass) was discharged by `task57`. | table-only (its one other inbound is `task57…:5`, itself retiring) |
| `task57_psi_bc_matvec_docs_pass.md` | RETIRE | Its own "distinctive lesson" is routed by the file to L-020, and digest §4 carries it verbatim: *"When the deletion is a COROLLARY of a design unification the SECTION'S THESIS is stale, not the line."* Optional micro-sharpening **S21**. | table-only |
| `task58_staleness_sweep_fix.md` | RETIRE | A per-file fix inventory whose one general lesson the file itself routes to the archive (*"LESSON REINFORCED → lessons L-021"*); its residual false-positive classes are digest §1d's `hasattr`-is-False-for-a-dataclass-field ladder, and its transient `:ref:`-to-a-title-less-anchor warning is digest §2a's `ref.ref` clause. | table-only |

### The two OPEN issues

`#240` and `#263` are open, which is the obvious objection to retiring
`feedback_d5b_d6_campaign_closeout.md` and `feedback_issue_257_s9_…md`. It does not hold. Both
memos are **docs close-outs for work that landed** — the D6 pass cleared all seven `.. todo::`
stubs, and the S9 criterion section is live at `operator_algebra.rst:5529`. What the open issues
track is the remaining *campaign*, and under Cardinal Rule 4 the issue is where that plan lives.
Neither memo is the plan; each is a report about a finished docs pass whose output is on the page.

---

## Salvaged

Each entry gives the lines to carry and the destination. **S3–S9 go to a new archive section
L-114** (next free number; the file has no `_archive/` directory — `lessons_archive.md` is the
cold layer). The rest are one- or two-line insertions into the named digest bullet, so the digest
grows by about a dozen lines against 4 200 retired.

*The numbering skips S11 and S15 deliberately: both started as salvages and became something
else during verification — S11 became refuted claim **F3**, S15 became skill uplift **U1**. The
gaps are left rather than renumbered so the table's references stay stable.*

### S1 — a cut can ORPHAN a citation, and that direction is gated
*From `feedback_phase2c_staleness_sweeps.md` §4; re-measured this session (M1), not relayed.*

> **Cutting a section can ORPHAN a citation, and unlike a dead code-xref that direction IS
> build-gated** — a `[Key]_` defined but no longer referenced warns
> `Citation [Key] is not referenced. [ref.citation]` at DEFAULT severity. Resolve by moving the
> citation to a surviving section covering the same topic, or deleting it when a sibling section
> already carries its own citations for that topic.

**Destination:** `lessons.md` §2a, appended to the "These DO warn, so they are gated" bullet. It
is the deletion-direction complement of §9's existing `grep '^\.\. \[Key\]'` clause, which covers
only the *adding* direction.

### S2 — a relocation's LoC target never subtracted the tables
*From `feedback_phase2d_table_preservation_constraint.md`.*

> **A relocation brief's headline LoC target is an estimate written WITHOUT subtracting the
> content it also tells you to preserve.** Count the preserved tables first
> (`awk 'NR==X,NR==Y' file | wc -l`); the achievable delta is `section − tables − stub`, and a
> shortfall against the headline is CORRECT, not a miss — state it with the clause that forced it
> rather than cutting further to hit the number.

**Destination:** `lessons.md` §1a, beside "A brief's named target can measure ZERO", as the
second face of the same rule (`process-discipline` "Measure a brief's premise before arguing its
scope").

### S3 — an operator's spelling sweep must chase the operator's *other* forms
*From `425_sn_chapter.md` §1.*

> **A spelling census is a FLOOR; the sites it cannot see decide the pass.** When a pass changes
> an OPERATOR, grep its **SPLITTING**, its **ITERATION MATRIX** and its **PRECONDITIONER**
> spelling too — `ψ_{n+1} = (L+C)^{-1}(Sψ + Bψ + q)` carries no `L+C−S−B` substring, so the
> census is blind, and leaving it makes the page state a five-member operator whose own splitting
> drops a term. Two more classes the census cannot see: a SECTION HEADING naming the count, and a
> non-rendered machine-facing `.. (vv-status rationale)` COMMENT restating the retired member
> list.

### S4 — a count-word names a SET, and the sets differ between adjacent pages
*From `425_sn_chapter.md` §2 and `425_outside_chapter.md` §3.*

> **When a count changes, write the MEMBERS beside it, or rename to a countless form.** One page
> called the algebra "five operators" meaning `{L,C,S,B,F}`, its neighbour "four operators"
> meaning `{L,C,S,B}` — two different sets, one numeral. And triage every count-word hit **by
> referent**: of 29 `(four|five|six)[- ](term|operator)` hits, most were a different four. The
> nastiest member is a number that matches nothing on the page (*"the four operator families"*
> over a matrix whose braces name five operators in three `BlockRole` groups) — **do not bump
> such a number; re-derive what it counts**, and if nothing counts it, drop it.

### S5 — scope a failed universal to the PASS, where the diff can check it
*From `425_sn_chapter.md` §3.*

> **When a chapter-wide universal you published turns out false, the honest replacement is a
> claim about the PASS, not about the chapter** — *"this pass changed no measured value; every
> edit is an algebra spelling"* is checkable from the diff, where *"every fixture in this chapter
> is Σ₂ₙ ≡ 0"* was false (an adjoint page carried a nonzero ladder and two gates inject one).
> Publish the real census beside it with its denominator and the reason the denominator is
> complete.

### S6 — the cross-solver claim nobody re-reads
*From `425_outside_chapter.md` §6. The strongest single item in the 28.*

> **A page that says "the other method does X instead" is the sentence nobody re-reads**: it
> lives on page A and its truth lives in solver B, so neither method's maintainer ever revisits
> it. One such sentence had been retired by an ERR entry two months before the pass and was still
> present tense. ⟹ when a page compares to another solver, **open that solver's function AND grep
> the error catalogue for the term** — the catalogue entry is where a reversal is recorded.
> Repair shape: a `.. note::` preserving the retired sentence verbatim in quotes, naming what
> retired it, and stating the surviving difference.

### S7 — the permissive-shim docutils parse mints findings the real build does not have
*The refutation of `425_outside_chapter.md` ⛔ item 2, measured this session.*

> **A permissive-shim docutils parse is a real instrument for structure (it catches
> `Unexpected indentation` and `Bullet list ends without a blank line` that no hand-rolled check
> sees) — but it FABRICATES markup findings for any directive the shim replaced with a
> nested-parsing stub.** `[M]` 2026-09-21: the claim that a `list-table` cell's continuation line
> beginning `-` "is a bullet and silently eats the role" does **not** reproduce under Sphinx —
> breaking a `:math:` before the operator renders identically to breaking it after, no warning,
> no `<ul>`. The only shape that does break is a continuation indented back to the **cell-marker**
> column, and that is LOUD, not silent: `ERROR: … uniform two-level bullet list expected …
> (3 vs 2)` plus a start-string-without-end-string warning. ⟹ the digest's "validate your own
> parser before believing its NEGATIVES" has a mirror — **reproduce a shim POSITIVE against the
> real builder before publishing or acting on it**. The memo half-knew this: it recorded six
> list-table anomalies that were all artefacts of its own crude table parser.

### S8 — a scope note is load-bearing only BESIDE the claim it scopes
*From `425_outside_chapter.md` §8.*

> **A correct explanation one paragraph away is read by nobody who lands on the equation.** Twice
> a page already carried the right scope note 2–6 lines below the site and it did not count. Move
> it up; that is a genuine improvement, not gate-gaming — and it is why an annotation window is
> specified in lines.

### S9 — the grey zone: date the spelling in place
*From `425_outside_chapter.md` §9; a sixth disposition for digest §4's five-register tense sort.*

> **Present-tense grammar with a historical subject is the genuine grey zone between "live
> guidance" and "history". Resolve it by DATING THE SPELLING IN PLACE** — *"the two-gain spelling
> here is Wave O's, which is what the cited captures were taken against; `N_{2n}` joined at CS4c
> step 3 and rides this argument unchanged"* — which keeps the evidence honest instead of
> retro-fitting a current member list onto a measurement that never saw it.
>
> Companion, from the same section: the (ii)-vs-(iii) test is **checkable**, not a judgement —
> read the composition site, or read the fixture's constructor call, where the **ABSENCE** of a
> kwarg (`sig_2` defaults to `None`) is the proof.

### S10 — a bare `Γ_-` in prose is an RST reference
*From `feedback_bc_trace_law_wave_12.md` pitfall 2; re-measured this session (M2).*

> **A symbol written in prose with a TRAILING underscore is RST reference syntax** — `Γ_-`, `S_-`
> and `X_+` all parse as a reference to target `Γ` / `S` / `X` and raise
> `ERROR: Unknown target name` at default severity; `-W` exits 1 on it. The control `V_bulk` (no
> trailing `_`) and the escaped `Γ\_-` are both silent. Fix by the math role in mathematical
> context or by a literal / escaped underscore in headings and prose — which `[M]` is already the
> live convention: all 7 corpus occurrences sit inside a literal, a code comment or a directive
> option.

**Destination:** `lessons.md` §2a, "These DO warn, so they are gated".

### S12 — the ARCHITECTURE pass has no SymPy, and its source ladder is different
*From `feedback_wave_o_operator_algebra_docs.md` §"Source-of-truth discipline".*

> **For an architecture pass there is no Branch-1 SymPy module, and the algebra of record is a
> different ladder:** the rich CODE docstrings (module / class / property) → the throwaway
> diagnostic instruments in `derivations/diagnostics/` (**READ and RUN them** to confirm every
> cited number) → the commit bodies (which carry the exact ULP figures and the rationale — QUOTE
> them) → the cross-domain-attacker frame memos, where the hardest "why" lives. Pull numerical
> bounds from the TEST FILES, never from a brief's memo estimate.

**Destination:** `lessons.md` §6, as the sibling of the existing L-005 read-order bullet (which is
the constructive-math ladder ending in SymPy). Also the natural home for a one-line note in
`algebra-of-record`, which currently addresses only the constructive-math case.

### S13 — the growing-family peer-section template
*From `feedback_wave_o_operator_algebra_docs.md` §"The shared section template".*

> **When a campaign lands step after step into ONE canonical page, every step gets a peer
> `====` section on a shared 8-part template:** (1) Lead — commit chain, issue, date, and one line
> placing the step relative to its predecessor; (2) Key Facts admonition, 5–7 bullets, always
> including one honest-scope bullet; (3) labelled equations, each with its `vv-status`;
> (4) ONE list-table that is the canonical at-a-glance index; (5) the load-bearing rationale,
> usually a rejected-alternative catalogue — this is the Cardinal-Rule-3 payload; (6) a
> numerical-evidence list-table; (7) an honest-scope `.. warning::` framed as the attacker's
> ABSENCE, not a hedge; (8) cross-references to the predecessor sections. The per-solver page gets
> only Key-Facts bullets pointing here — never a duplicated derivation.

**Destination:** `lessons.md` §6, one line with the eight parts; it is the "growing family"
sibling of AGENT.md's 9-step close-out arc and of the live
`feedback_canonical_convention_page.md`'s one-page anatomy.

### S14 — tombstone the SCAN PATH, not only the deep site
*From `feedback_err058_success_closeout_supersedes_phase_chain.md`.*

> **A deeply-nested phase chain has two reading paths and both need a tombstone.** A reader
> scanning the phase block's top-of-section step summary hits the stale terminal claim ~700 lines
> before the deep tombstone ever renders. ⟹ one forward-pointer `.. note::` at the TOP of the
> phase block listing which terminal decisions were reverted, **plus** inline tombstones at each
> deep site. A single close-out section does not catch a scanning reader.

**Destination:** `lessons.md` §4, appended to "Preserve the WHY; tombstone, don't delete".

### S16 — only a too-SHORT underline warns
*From `feedback_d5b_d6_campaign_closeout.md`; re-measured this session (M3).*

> **An over-long section underline is silent; only `Title underline too short` warns.** So the
> underline scan has one direction that matters, and normalising an over-run is cosmetic — scope
> it to your own lines and never "fix" a pre-existing one.

**Destination:** `lessons.md` §9, appended to the length-changing-rename underline bullet.

### S17 — the marker-ladder error points at a section you did not touch
*From `feedback_issue_247_legA_mode10_closeout.md`; the diagnostic half AGENT.md lacks. Its
positive control is `feedback_issue_251_legB…`, where reusing the file's marker gave zero
collisions first try.*

> **RST marker levels are file-local AND assigned by FIRST-APPEARANCE ORDER, so introducing a
> never-before-seen marker EARLIER in the file than an existing same-depth marker demotes the
> existing one** — a new `"""` at line 3267 claimed level 5, pushing the file's existing `'''` at
> line 5865 to level 7. ⚠ **The resulting `Inconsistent title style: skip from level 5 to 7`
> ERROR points at the OLD sections, not at your edit** — do not read it as "the old section is
> broken". Before adding a level the file has never used, grep the file for an existing marker at
> that depth and REUSE it.

**Destination:** `lessons.md` §9, appended to the marker-ladder bullet; AGENT.md's
"Title markers, labels, and vv-status" carries the rule and should gain the diagnostic clause.

### S18 — `| tee` reports tee's exit, not sphinx's
*From `feedback_phase_c_stub_expansion.md`.*

> `python -m sphinx -W … | tee log` returns **tee's** exit status, so a failing build reads as a
> pass. Redirect instead (`2>/tmp/log; echo $?`).

**Destination:** `lessons.md` §2c, folded into the existing
"⚠ `nohup … &` inside a background Bash call reports the SHELL's exit, not sphinx's" bullet — same
failure, second member, so the bullet becomes about wrapper-exit generally.

### S19 — name the FORBIDDEN SENTENCE
*From `feedback_issue_257_s9_coherent_promise_criterion_seam.md` LESSON 3. Squarely in the
Directive-5 curator lane.*

> **When the whole point of a section is that the obvious reading is WRONG, add an explicit
> `.. warning::` naming the forbidden sentence verbatim and giving the correct framing.** A
> future session quoting the page for V&V reasoning meets the warning before the misreading. It
> is the page-level analogue of `vv-principles`' *"NEVER write 'MMS verifies the eigenvalue'"*,
> and it is stronger than a hedge because it is greppable.

**Destination:** `lessons.md` §7, beside the existing coverage-claim bullet.

### S20 — the page as COUNTER-RECORD to a surviving wrong artefact
*From `feedback_phase_d_carlson_seed_narrative.md`.*

> **When a fix OVERRIDES an upstream artefact that still exists and will be read first (a
> literature memo's implementation note, a plan's pseudocode), the page is not merely recording
> the outcome — it is the permanent COUNTER-RECORD.** Give the correction its own named
> subsection and carry the falsifying evidence table verbatim (here: four interventions ×
> site × residual, with the no-op, the pass, the redundant and the degenerate-coincidence rows).
> Without the table the next session re-applies the wrong injection point, because the surviving
> memo says so. This is the doc-side of `process-discipline`'s "a refuted candidate is
> first-class output": the structural reason must be findable *before* the artefact it refutes.

**Destination:** `lessons.md` §1f, beside "Keep a refuted prediction — the refuted MECHANISM is
the interesting half".

### S21 — optional micro-sharpening
*From `task57_psi_bc_matvec_docs_pass.md`; sharpens digest §6's "A correction sweep must not
acquire a SECOND SUBJECT" by naming the boundary.*

> Co-fix a neighbouring defect **only where it sits INSIDE a clause you are already rewriting**;
> a standalone instance of the same family is flagged for its own pass, with its sites listed.

---

## Referrer edits

Five, all specified in full. Every other inbound line lives inside a retiring candidate and
disappears with it — stated explicitly per row in the table above so nothing is assumed.

**R1 — `.claude/agent-memory/archivist/feedback_canonical_convention_page.md:141-142`**
(a surviving file; re-point **and** correct, since the cited rule has since been sharpened)

- current: ``The acceptance gate is **baseline-warnings-unchanged** (per the`` /
  `` `feedback_bc_trace_law_wave_12` rule), NOT count=0. ``
- replace with: ``The acceptance gate is **baseline-unchanged** (`AGENT.md` § "The build gate"), ``
  ``NOT count=0 — and the gate is the WARNING/ERROR/CRITICAL **set**, freshly measured with `-E` ``
  ``each session, never a count and never a quoted baseline. ``

**R2 — `.claude/agent-memory/archivist/feedback_capstone_architecture_page.md:11`**
(a surviving file; re-point)

- current: `[[feedback_post_wave_cleanup_docs]] (close-out arc) but distinct: this`
- replace with: `the post-wave close-out arc (AGENT.md "Close-Out Narrative Arc"; its`
  `capability-flip variant is lessons.md §6 → L-076) but distinct: this`

**R3 — `.claude/agent-memory/archivist/feedback_curvilinear_aniso_norm_reconciliation.md:136`**
(a surviving file; rewrite as a pointer to the durable record)

- current: `Cross-ref: [[feedback_err058_success_closeout_supersedes_phase_chain]]` /
  `(the #195 close-out this builds on), …`
- replace the wikilink with: `the ERR-058 / #195 close-out this builds on (#195 CLOSED; the`
  `record is `docs/theory/verification/error_catalog.rst` ERR-058 plus the SN curvilinear pages)`

**R4 — `.claude/agent-memory/archivist/feedback_curvilinear_aniso_norm_reconciliation.md:137`**
(same file, next line; rewrite as a pointer)

- current: `[[feedback_issue_196_eigenvalue_verification_closeout]]` / `(the #196 sequel), …`
- replace with: ``the #196 sequel (#196 CLOSED; the record is``
  ``:ref:`sn-issue-196-bit-identical-vs-floor` in``
  `` `docs/theory/methods/sn/curvilinear_numerics.rst`) ``

**R5 — `.claude/agent-memory/method-implementer/issue_168_phase_d_closeout.md:229` — CROSS-OWNER**

I do not own this file and did not edit it. The line is a source-list bullet:
`- **Archivist feedback**: `.claude/agent-memory/archivist/feedback_phase_d_carlson_seed_narrative.md``.
Proposed edit for the method-implementer owner (or the orchestrator, applying it in the same
commit): **delete the bullet**, since the Phase-D narrative it points at is on
`docs/theory/methods/sn/` and in `error_catalog.rst` ERR-026/ERR-058, and the pass's transferable
lesson moves to `lessons.md` §1f (S20). Re-pointing it at the archive section would also be
correct; deletion is cleaner, because the neighbouring bullets in that list already name the
Sphinx narrative and the catalog entry directly.

---

## Refuted claims carried by the candidates

Each is a reason to retire rather than merely an absence of value: the file states, in the present
tense, something the corpus now knows to be false.

- **F1 — `feedback_phase2b_label_dense_partition.md`**, counter-pattern: *"Anchors used in `:eq:`
  form … won't be caught by Sphinx `-W`."* Refuted by digest §2a, which records a dangling `:eq:`
  as measured-with-controls WARNING, and by L-113 for `:ref:` / `:doc:` cross-doc. Only the
  `verifies`-marker half survives, and that half is already digest §3 plus AGENT.md's INFO-severity
  note.
- **F2 — `425_outside_chapter.md`** ⛔ item 2: the `list-table` continuation-line bullet trap.
  Refuted this session with two-sided controls; the refutation is salvage **S7**.
- **F3 — `feedback_bc_trace_law_wave_12.md`** pitfall 1: *"these are ERRORs not WARNINGs but
  `sphinx-build -W` does NOT promote them to failure."* `[M]` a docutils `ERROR: Unknown target
  name` gives `sphinx -W` **exit 1**, against a clean-control page's exit 0.
- **F4 — `feedback_phase_d_carlson_seed_narrative.md`** §"Equation-label discipline": *"the labels
  are identical strings, so the cross-document graph walks resolve to a single equation node
  either way"*, offered as the reason there was no duplicate-label collision. `[M]`
  `orpheus/sn/sweep/psi_half_angle_seed.py` carries **3** docstring `:label:`s and is
  `automodule`'d **nowhere** in `docs/`, so the real reason is that the module is never rendered.
  The rule as stated breaks the moment anyone surfaces the module — which is exactly AGENT.md's
  standing warning that a module with `:label:` docstrings must be cross-referenced in prose
  instead.

Two further files carry an instruction that violates a standing rule rather than a false claim:
`feedback_err058_…` prescribes `git stash` for baselining (forbidden by `process-discipline`, and
replaced by digest §9's `git archive HEAD` into a temp tree), and `feedback_phase2c_…` /
`feedback_phase2b_…` prescribe a per-step or per-group build cadence that digest §9's
build-TWICE sequencing supersedes.

---

## Skill uplifts proposed (AGENT.md Directive 5)

Neither is a memory salvage; both are edits to a skill, raised here for the orchestrator to route.

**U1 — `vv-principles`, § "Bit-identity vs principled-equivalence": the ENTRY POINT decides which
claim class is even available.**

The skill gives three conditions for accepting a non-bit-exact change but never says which claim
is reachable. Proposed addition:

> **The entry point decides the claim class, so ask it before writing "bit-identical".** A
> fixed-source solve (one operator, one quadrature, `L.solve` vs Krylov-on-`apply` = the same
> `L⁻¹` arithmetic) CAN be bit-identical. The same two inner solvers under an eigenvalue entry
> point are wrapped in power iteration and are two different iteration schemes converging to the
> same fixed point only to the inner tolerance: **floor-equivalent, the same physics, not the
> same arithmetic** — never call it bit-identical. **tell:** a blanket "SI ≡ Krylov bit-identical"
> whose fixture list spans both entry points.

`[M]` all three entry points are live (`orpheus/sn/solver.py:2264` `solve_sn`, `:3241`
`solve_sn_fixed_source`, `orpheus/numerics/eigenvalue.py:374` `power_iteration`), and the corpus
already states the distinction at `docs/theory/methods/sn/curvilinear_numerics.rst:2995` — so this
uplift moves a shipped, anchored ruling into the vocabulary that `qa`, `test-architect` and
`numerics-investigator` read, which is the gap Directive 5 names.

**U2 — `algebra-of-record` (or `vv-principles` anti-patterns): analogical generalisation of a
closure formula.**

The one proposal of `peierls_greens_phase123_rich_narrative.md` that has **not** landed (`[M]` no
skill file matches `analogis|analogiz|parameter substitution`). Its catalogued instance is
ERR-035 (`error_catalog.rst:2433`, with catchers). Proposed wording:

> **NEVER analogise a closure formula across geometries by parameter substitution** — derive the
> new geometry's closure from first principles via Branch 1 on a *non-uniform* source, and
> cross-check structurally-independently against the rank-N first-principles form. **tell:** a
> closure for geometry B obtained by substituting B's chord length into A's formula.

Its two siblings are already covered and need nothing: uniform-source-test blindness by
`vv-principles` Mode 7 (which requires every test to declare which terms its ansatz activates and
which it nulls), and convergence-rate-as-bug-fingerprint by `numerical-bug-signatures`
Signatures 6 and 7, in sharper discriminator form than the memo proposed.

---

## NEEDS:

1. **Nexus was not consulted.** The brief asked for `context("<node id>")` on one candidate if the
   tools were available. They surface in my tool list, but I did not use them for the referrer
   census: the census is the load-bearing claim of this audit, `nexus-tools` warns that a graph is
   a build-time snapshot, and `git status --porcelain -- docs/` being clean does not tell me when
   the graph was last built. I ran two independent text censuses instead (the brief's and my own,
   different predicates, agreeing exactly, with a positive control). If the orchestrator wants the
   graph-edge view as a third instrument, `context()` on
   `.claude/agent-memory/archivist/feedback_bc_trace_law_wave_12.md` is the one candidate with
   three inbound edges and therefore the best test of whether the graph indexes agent memory at
   all.
2. **R5 is cross-owner and I did not apply it.** `.claude/agent-memory/method-implementer/` is the
   method-implementer's; the edit above is a proposal, and if that owner is being audited in the
   same campaign it should land in their commit rather than mine.
3. **L-114 is a number, not a reservation.** `lessons_archive.md` currently ends at L-113. If
   another owner's proposal also claims L-114 in this campaign, the archivist section should take
   whichever number is free when the commits are ordered.
4. **Not obtained: why #425's pass never got an archive section.** I established the gap (`[M]`
   L-099 → L-103 skips 2026-09-07) but not its cause. It may be that the pass ran with two
   concurrent archivist instances and neither owned the write-up. Worth one line in the campaign
   close-out if the orchestrator knows, because the same dispatch shape will recur.

---

**APPLIED 2026-09-22, by the archivist, at the orchestrator's instruction** (proposal accepted in
full after review with session context; the orchestrator commits). All 21 salvages landed in the
`lessons.md` bullets named above; S3–S9 landed as `lessons_archive.md` § **L-114** (230 lines,
`MEMORY.md`'s next-free pointer bumped to L-115); S17's diagnostic clause landed in
`.claude/agents/archivist/AGENT.md` § "Title markers, labels, and vv-status"; referrer edits
R1–R5 applied, R5 included (the cross-owner one-bullet deletion in
`.claude/agent-memory/method-implementer/issue_168_phase_d_closeout.md`); the 28 candidates are
`git rm`-ed; U1 landed in the SOURCE `docs/development/skills/vv-principles.md` with
`budget_tokens` 12700 → 13000 (measured ≈12990) and `tools.harness --check` reporting
`24 targets, 0 problems, 0 drifted; citations by ID: 726 resolved, 0 dangling`; U2 landed in the
hand-maintained `.claude/skills/algebra-of-record/SKILL.md` § "For method-implementers", citing
ERR-035. ⚠ **One referrer BOTH censuses missed** was found by `git grep` after the `git rm` and
repaired as **R6**: `.claude/scratch/open_fronts_audit.md:78` cited
`feedback_phase2c_staleness_sweeps.md` by bare filename. The mechanism is worth keeping — both
censuses excluded a path matching `scratch/`, and `.claude/scratch/` is a TRACKED directory that
matches it, so the two instruments were not independent on that axis and agreed tautologically
about what the shared exclusion hid (`instrument-doctrine` X4; the exclusion is a predicate, X2).
`docs/theory/conventions/indexing_and_layout.rst:663,2186` also matched a loose grep but is a true
negative: it names the **method-implementer's** `issue_196_pr_index_*_closeout.md` glob, 6 of whose
files exist. No Sphinx build was run and `docs/theory/verification/error_catalog.rst` was not
touched, per instruction.
