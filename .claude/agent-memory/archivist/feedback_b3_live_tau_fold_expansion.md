---
name: b3-live-tau-fold-expansion
description: "#236 Phase 2 B3 — expand the sn-tau-c-on-cellvisit-live stub into rich narrative. 4th in the A/B1/B2/B3 carve chain on discrete_ordinates.rst. A bit-identical (0-ULP) carve moving the SOURCE of τ from geometry-StreamingTerms.tau_mm to the angular closure; the precondition for Step C (geometry-τ retirement)."
metadata:
  type: feedback
---

#236 Phase 2 B3 stub→rich expansion in
`docs/theory/discrete_ordinates.rst` (`.. _sn-tau-c-on-cellvisit-live:`).
4th of a 4-section carve chain (A=`sn-tau-closure-owned`,
B1=`sn-closure-c-constants-owned`, B2=`sn-closure-c-on-cellvisit`, then
this B3). The siblings A/B1/B2 are STILL `.. todo::`-wrapped rich stubs;
B3 is the one elevated to full production-depth narrative. A
bit-identical CARVE (moves the SOURCE of an already-correct number, not
the number) — like affine-in-σ, the narrative's job is the WHY of a
re-sourcing, not a bug.

**Why (source order, algebra-of-record):** the method-implementer
CLOSEOUT MEMO
(`.claude/agent-memory/method-implementer/issue_236_phase2_b3_tau_c_live_closeout.md`)
was the load-bearing source (6-seam map, the fixture-ripple lesson, the
in-process 0-ULP proof, the Mode-11 τ-arm). The crosswalk plan §9
(`issue_236_phase2_tau_carve_crosswalk.md`) gave the L20 6th-consumer
finding (sweep_cache→loss_representation.py:3270 scan split) the §8
5-seam framing missed. Then READ the 6 production seams to learn the
actual landed shape (docstrings on `CellVisit.tau`, `_make_cell_visit`,
`cell_balance_terms`, `DD.update`, `sweep_cache.from_mesh_and_quad`,
`pole_angular_closure._build_per_ordinate_cache`/`tau_per_ordinate`) —
the docstrings ARE the prose seed.

**The 6-H4 shape that worked (H3 title KEPT, sub-sections at `^^^^`):**

1. **Lead framing** (no header) — situate in the A/B1/B2/B3 chain (B1=cache
   populator, B2=residual twin carries c, B3=LIVE paths); state
   precondition-for-Step-C up front; 3rd-of-4-c-sites + 5th-τ-consumer;
   bit-identical = moves SOURCE not number; roadmap the 6 sub-sections.
2. **"τ is an angular-scheme property the closure owns"** — the WHY:
   `mm-weights` is built ENTIRELY from (μ_n, μ_{n±1/2}, level partition),
   NONE of which is spatial-streaming geometry. τ attaches to an
   ORDINATE not a CELL. Historical accident (factory was in scope).
   Closure is the owner; `tau_per_ordinate` on Protocol AND ABC (B1
   both-sites mint); MM vs Identity dispatch BY TYPE no `coord==`;
   `_build_c_per_ordinate_cache`→`_build_per_ordinate_cache` (3 constants).
3. **"The spatial ⊗ angular separation forbids coupling DD to the
   closure"** — the architectural crux: DD is STATELESS, reads only
   CellVisit+UpstreamState. Coupling DD→closure collapses the tensor
   product into a Cartesian-vs-curvilinear conditional = the dispatch the
   unified body deleted. So constants travel as DATA on CellVisit
   (c_in/c_out B2 + tau B3), stamped at ONE site `_make_cell_visit` (all
   4 dag_walk paths funnel, Pattern 2). Include the actual stamp
   code-block; spatial axis NEVER reaches the angular axis via the type
   system.
4. **"Why the CellVisit.tau default is 1.0, not 0.0"** — τ=1 is the
   NEUTRAL element: recurrence `(ψ̄−(1−τ)ψ_in)/τ →[τ=1] ψ̄` (identity,
   exactly slab's no-redistribution need); c_out, 1/τ all well-defined.
   0.0 = ÷0 LANDMINE (every consumer divides by τ). 1.0 default makes an
   un-stamped visit compute the IDENTITY (safe no-op) → a missing stamp
   surfaces as wrong-but-FINITE (snapshot catches) not nan-propagation.
5. **"The three live consumers and the L21 framing"** — (1) cell_balance_terms
   solve dir (now reads c only, τ NOT AT ALL — right factoring: denom
   needs (ΔA/w)c_out, numer needs (ΔA/w)c_in·ψ); (2) DD.update angular
   recurrence (MINT eq-label `dd-mm-angular-recurrence`, reads visit.tau,
   the line the 1.0 default protects); (3) CumprodScan split (MINT
   `dd-mm-scan-split` tau_inv/mm_a_in_coeff, consumed at
   loss_representation.py:3270, algebraically ≡ recurrence (2) — L21 twin,
   two applications of ONE operator). Pattern 5: closure exposes only
   PRIMITIVE τ; consumers derive 1/τ, (1−τ)/τ, α_out/τ locally.
6. **"Why the fold is bit-identical, and the regression floor for Step
   C"** — 2-leg argument (closure-τ 0-ULP=geometry-τ via Leg-1 gate +
   per-level→(N,) gather is pure permutation). In-process |visit.tau −
   st.tau_mm|≡0 at sphere(single-level)+cyl(multi-level). Then the THREE
   committed catchers = regression floor for C: Leg-1 (values),
   test_cell_visit_c_stamp τ-arm (ordinate MAP, vv L11 Mode-11 — the
   ×1.1 mutation drifts cyl scalar 0.2% w/ no OTHER red, structurally-
   independent geometry-vs-closure producers), test_affine_carve_baseline
   (seam-6 scan). Step C deletes spherical/cylindrical_streaming +
   StreamingTerms.tau_mm/alpha_in/out + tau_mm_per_level confident
   nothing live reads them.

**Disciplines confirmed (this page):**

- **No new CLAIM :label:.** The carve sharpens existing claims
  (`mm-weights` count 75, `dd-solve` 7, `dd-mm-closure-constants` 3 — all
  verifies-targets, CITE with `:eq:` NEVER re-`:label:`). The 2 NEW labels
  (`dd-mm-angular-recurrence`, `dd-mm-scan-split`) are DERIVATION/
  representational (the algebra OF the M-M recurrence) — added UNTAGGED
  (no `.. vv-status:`), sitting inside the verified narrative chain.
  matrix.rst AUTO-REGENs them into its orphan list (645-646), NOT a
  hand-edit, NOT a warning. Same posture as affine-in-σ's 4 derivation
  labels.
- **L21 cross-ref: pick the SN-context anchor, not the LD one.** First
  draft cited `:ref:`ld-ubld-unified-matvec`` (multi-D LD moment context)
  — but (a) that anchor doesn't exist (real one is
  `ld-ubld-unified-moment-matvec`) AND (b) wrong context. The right SN
  apply/sweep-twin anchor is `phase-c-sweep-frame-matvec` (the curvilinear
  apply-direction matvec, the genuine L21 twin of the curvilinear sweep).
  Intra-doc dangling `:ref:` IS caught by `-W` → the typo would have
  failed the build; grep `^\.\. _<target>:` BEFORE citing.
- **B3 stub was ALREADY rich (not a bare todo).** The method-implementer's
  stub had a full descriptive narrative INSIDE the `.. todo::`. The
  expansion's value-add over the stub: the WHY-derivations (τ-is-angular
  / ⊗-separation / 1.0-default), the L21 framing tied to the existing
  page anchor, the math eq-labels for recurrence (2) and split (3), the
  regression-floor framing of the 3 catchers, the stamp code-block. Don't
  just un-wrap the todo — ELEVATE to production-section depth.
- **Build gate**: `-E -W --keep-going` baseline = 1 (the standing
  `Mesh1D.from_geometry :paramref:` ERROR). Post-edit identical = 1.
  Count-unchanged, NOT zero. MAIN-checkout baseline is 1 (not the
  worktree's 11).
- **File ladder = `=`/`-`/`~`/`^`.** H3 (the B3 carve sections) is `~~~~`;
  sub-sections inside go to H4 `^^^^` (file HAS `^^^^` at 1974/1997/2054).
  KEEP the existing H3 title + its 63-byte/59-codepoint `~` underline
  (em-dash 3 bytes, τ 2 bytes; ≥ title len fine).

**Quality self-assessment (1-5):**

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | full WHY-derivations: τ-is-angular from mm-weights inputs, ⊗-separation, τ=1 neutral-element reduction, the recurrence≡scan-split algebra |
| Cross-references | 4 | every :eq:/:ref:/:func:/:attr:/:meth:/:mod: resolves; :class:/:meth: render plain-text by page convention (not my regression) |
| Numerical evidence | 5 | in-process \|visit.tau−st.tau_mm\|≡0 sphere+cyl, 588+60 snapshots, the ×1.1 mutation 0.2% drift, 0-ULP 2-leg argument |
| Failed approaches | 5 | the 0.0-default ÷0 landmine, the coupling-DD-to-closure anti-move, the Mode-11 stamp-catcher gap (named twins never call the rewired reader) |
| Code traceability | 5 | all 6 seams cited by symbol, the stamp code-block, the 3 catchers named, Step-C retirement targets enumerated |
| Derivation source | 5 | closeout memo + crosswalk §9 + the 6 production seam docstrings (the algebra-of-record for a re-sourcing carve is the operator code) |

Weakest = cross-references (4), only because the page-wide plain-text
code-xref convention (`orpheus.sn.spatial.*` not member-autodoc'd with
resolvable anchors) is out of scope for a stub expansion.
