---
name: phase-d-carlson-seed-narrative
description: Phase D Carlson coupled-pole seed expansion (Issue #168 Phase D) — pattern for documenting a 4-memo-input narrative that closes ERR-026's identity-and-rate sub-claims while narrowing the open scope. ~880 LoC added across discrete_ordinates.rst (primary) + boundary_conditions.rst (companion). The literature memo's architectural sketch was empirically falsified by the diagnostic memo; the docs MUST preserve the corrected-injection-point story as load-bearing pedagogy.
metadata:
  type: feedback
---

When a session ships a fix that **overrides** the original plan's
architecture based on empirical diagnostic evidence (literature memo
proposed site X; diagnostic memo proved site X was a no-op and site Y
was canonical), the rich narrative MUST preserve the corrected-
injection-point story as a named subsection. Phase D shipped exactly
this: the literature memo routed Carlson seed into WDD ``psi_face_in``;
the diagnostic falsified that with intervention `[A]` and confirmed
intervention `[B]` (M-M angular recurrence's ``psi_half_left``).

**Why:** Future sessions reading the closeout (in 6 months, 1 year)
will re-read the literature memo's §7 implementation note FIRST and
re-apply the wrong injection point. The Sphinx narrative is the
**permanent counter-record** — the diagnostic table with all 4
interventions (`[A]` FAIL no-op, `[B]` PASS, `[C]` redundant, `[D]`
degenerate-coincidence) IS the canonical evidence chain.

**How to apply:**

1. **Cite the diagnostic memo** in a named subsection ("The corrected
   injection-point story" / "X plan correction"). The user marked the
   diagnostic memo's 4 plan corrections as load-bearing in the brief
   — they expect them surfaced in Sphinx, not just in agent memory.
2. **Include the empirical 4-row table** (intervention × site ×
   residual) verbatim. The pattern (`[A]` no-op + `[B]` PASS + `[D]`
   coincidence-on-degenerate-probe) is the diagnostic mark — without
   it the next session will re-attempt `[A]` because the literature
   says so.
3. **Add the structural-independence cross-check** (vacuum-BC probe
   producing a non-trivial `phi_aux` profile distinct from `ψ_cell`)
   as a separate paragraph below the intervention table. This proves
   `[B]` is not coincidentally equal to `[D]` on a less degenerate
   probe.
4. **Name the falsification test explicitly** —
   `test_l1_carlson_vs_zero_seed_differ_on_vacuum_bc_probe` — and
   explain why a flat-ψ-reflective probe alone cannot catch a future
   regression that replaces Carlson with broadcast-cell-centre.

### Equation-label discipline: docstring vs RST source-of-truth

The method-implementer declared 3 `:label:` blocks in the
`psi_half_angle_seed.py` module docstring (`hebert-3-434`,
`hebert-3-435`, `hebert-3-432-source`). The brief instructed:
*"Add the corresponding `.. math:: :label:` blocks in the theory
page proper so they exist as graph-walkable equation nodes."* This
creates a potential duplicate-label collision because Sphinx renders
both the docstring and the RST page.

**Resolution:** Make the RST page the **canonical declaration**, and
the docstring `:label:` blocks become "module-local declarations
mirroring the canonical RST source-of-truth". Since the labels are
identical strings, the cross-document graph walks resolve to a
single equation node either way. The `.. note::` block in the RST
narrative explicitly states this — "the module docstring is the
canonical algebra-of-record; the Sphinx page is the presentation
layer".

**Why:** Per `algebra-of-record` skill, the SymPy / Python derivation
module is the source-of-truth, the Sphinx page narrates it. A
duplicate-label collision would force one to be hidden; the
narrative resolution avoids the collision because the labels are
**identical** (same string, same equation, same canonical anchor).

### Phase D's ERR-026 PARTIAL → PARTIAL narrowing

The brief explicitly warned against "Phase D is complete" framing
and required PARTIAL → PARTIAL with three-sub-claim breakdown
(identity CLOSED, rate CLOSED, magnitude OPEN per Issue #195). The
canonical structural pattern:

* **3-row sub-claim table** with status + closed-by columns. Each
  row maps to a Wave's contribution. This is the load-bearing
  evidence for the `error_catalog.md` ERR-026 status update.
* **Pre-asymptotic-magnitude paragraph** explaining WHY the rate is
  correct but the magnitude is not yet ≤ tolerance at practical nx.
  The candidate follow-up paths (higher-order pole-face spatial
  closure, or `phi_0` recomputation via M-M output) are the
  open-research-paths content per the 9-step CLOSED arc adapted to
  the PARTIAL variant.
* **Xfail-strict markers stay** — they will xpass under the Phase D
  defaults, triggering the deferred Step 5 marker-removal commit;
  but the pre-asymptotic-magnitude regression prevents flipping to
  `strict=True` until Issue #195 lands.

This is the **partial-OPEN close-out** variant of the 9-step arc
(AGENT.md §"Close-Out Narrative Arc", the partial-OPEN variant). The
ERR-026 entry stays PARTIAL through Phase D Step 3; only Step 5
closes the markers, and Issue #195 closes the residual magnitude
question.

### Architectural-choice 2x2 framing

When a design decision rejects an alternative (Option α composition
vs Option B sibling Protocol), the Sphinx narrative MUST include
both options with rejection rationale. Phase D's "Option α (shipped)
vs Option B (rejected, SRP violation)" subsection is the canonical
form: 2 bullets, each with 2-3 lines of why-shipped or
why-rejected. This is for the future session who reads the
production code and asks "why isn't this a separate Protocol on
SNMesh?" — the rejection rationale is in the docs, not in
git-blame archeology.

**The L = 0 isotropic-only WARNING block** is the
load-bearing future-refactor safeguard. The class docstring already
carries a WARNING admonition (per the closeout memo's deviation #2);
the Sphinx narrative re-states it in the same WARNING admonition
form. Mode 6 (convention drift) is the failure mode flagged.

### Companion BC section pattern (boundary_conditions.rst)

When Phase D's `bc.apply` call sequence changes from once-per-matvec
to twice-per-matvec, the BC trace contract page needs a **table
showing the two calls** (caller / input / output use). The pattern:

* Update the existing §16A.3 admonition note to say "at least
  once per matvec" + cross-link to the new Phase D section.
* New section `bc-phase-d-two-bc-applies-per-matvec` with the 2-row
  call table + the capture-and-compare strengthening explanation.
* Locate-by-shape-AND-content discipline in the test: the gate's
  matching strategy must protect against a future third BC apply
  with the right shape but different semantics.

The two-application-sequence is structurally distinct: Call #1 is
a Phase D linear-in-ψ extraction (formally not a §16A.3 trace input
since cell-centre values flow in, but the resulting scalar is
linear and matches the canonical inward-zero-weight ordinate's
value on the load-bearing test configurations); Call #2 is the
canonical Phase C trace law. The narrative MUST be explicit that
Call #1 is a "principled shortcut" — not a contract violation, but
a linearly-compatible scalar extraction whose values match the
canonical form on production-relevant inputs.

### Build-warning gating discipline

Pre-existing baseline (Phase C era): 7 WARNINGs + 9 ERRORs = 16
total issues. The Phase C inconsistent-title-style errors at lines
2423-2656 are **pre-existing baseline noise** — they were shipped
with Phase C and stay through Phase D edits. The acceptance gate is
**count unchanged from pre-edit**, NOT count = 0.

The one new issue Phase D temporarily introduced was a duplicate
`[Hebert2009]` citation definition (the citation already existed in
`collision_probability.rst`; Sphinx resolves citations cross-document
but adding a NEW definition produces a `duplicate citation` warning).
**Resolution:** delete the new definition; the existing one resolves
through Sphinx's cross-document citation graph. This is the same
pattern as [[feedback_phase2c_staleness_sweeps]]'s "verify
[Sanchez1986] etc don't exist before adding citation marker" —
grep BEFORE adding any `.. [CitationKey]` block.

### Quality score self-assessment

| Dimension              | Score | Justification                                    |
| ---------------------- | :---: | ------------------------------------------------ |
| Derivation depth       |   5   | Hébert (3.432)–(3.435) full derivation; flat-ψ algebraic verification trace; vacuum-BC trace |
| Cross-references       |   5   | Every :func:/:class:/:mod: linked; cross-doc :ref: into boundary_conditions; back-refs to Phase C anchors |
| Numerical evidence     |   5   | 4-row intervention table; 12-cell Gate 1.1 crosstab; vacuum-BC profile array; convergence rate [3.33, 2.46] |
| Failed approaches      |   5   | Intervention `[A]` falsification + `[D]` coincidence; Phase B docstring justification falsified with mechanism |
| Code traceability      |   5   | Every equation tied to specific module/file/line; closeout memo cited; diagnostic script named |
| Derivation source      |   5   | psi_half_angle_seed.py module docstring is canonical algebra-of-record; RST narrates it |

Total: **30/30**. The ~880 LoC sweep is maximum-effort per Cardinal
Rule 3.
