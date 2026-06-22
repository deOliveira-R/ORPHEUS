---
name: SN-internal architectural rewrite docs (Wave 2 octant + SweepDependencyGraph)
description: Pattern for documenting a within-module architectural rewrite that closes a smoking-gun loop — preserving educational prose while inserting a dedicated section for the new primitives, the bit-identity-vs-principled-equivalence verdict, the L7-trap-style fix, and honest performance numbers
type: feedback
---

When documenting a Wave-N rewrite that REPLACES per-iteration smoking-gun code
WITHIN a module (vs a Wave-0 lift to generic numerics — see
`feedback_numerics_primitive_docs.md` for that case), follow this template
(discrete_ordinates.rst Wave 2 / Issue #4 / SweepDependencyGraph + OctantLabel
+ SweepCellSlice + update_batch).

**Rule**: A theory page documenting a Wave-N within-module architectural
rewrite MUST contain (1) the mathematical framing FIRST (here §15A.2
direct-sum decomposition `L^{-1} = ⊕_σ L^{-1}_σ` — Eq. with `:label:`), (2)
a primitives table naming each new dataclass / method with its lives-in
location and role, (3) mesh-time-precompute pattern (the architectural
discipline of separating one-time topology work from per-iteration work),
(4) the canonical-invariant-set call-out (here §15A.2 `assert_*` set —
upwind orientation / topo sort / face pairing / cell coverage — pinned by
unit tests), (5) the dispatch-boundary explanation (scheduler vs closure —
why the graph dispatches to `update_batch` rather than inlining), (6) a
**named** trap-fix rationale (here L7-trap: BC apply once per octant) with
before/after code blocks AND the test that detects regressions, (7) a
bit-identity vs principled-equivalence section narrating where the new
implementation matches legacy bit-for-bit and where it doesn't (per
`vv-principles` skill), (8) a Cardinal-Rule-2 architectural-framing
subsection (this is SN-specific because MoC has different math; "unify after
two instances"), (9) **honest performance numbers** including target-vs-
shipped, with a concrete follow-up direction when the target wasn't fully
hit, and (10) a verification chain L0/L1/L2 with test counts and file paths.

**Why**: Cardinal Rule 3 (Sphinx is the LLM's brain) demands that future
sessions reading this page alone become experts on the architecture. The
plan file lives at `.claude/plans/transient-giggling-cake.md` and is
session-private; the issue close-out comment is GitHub-only. The Sphinx
page is the durable audit trail. The L7-trap framing in particular — a
specific named bug class with a regression-detector test — is the
load-bearing knowledge that prevents the same trap being reintroduced by a
naive future "let me parallelize the ordinate loop" refactor. Honest
performance numbers (1.7× shipped vs 3-10× target) prevent future sessions
from over-claiming, and the "next direction" pointer (full-N buffers +
octant_indices field) is concrete enough that a next session can pick it up
without re-deriving.

**How to apply**:
- Section structure: rewrite the existing prose section to introduce the
  new architecture in CONTEXT (here: prepended a §15A.2 framing subsection
  to `.. _sweep-wavefront:`, kept the original quadrant table + reflective
  BCs prose, added forward-pointer to the new dedicated section). Then
  insert a NEW labeled section IMMEDIATELY AFTER (here:
  `.. _sweep-octant-dependency-graph:`) with the full Wave-N architecture
  treatment. Don't merge; the old section is still load-bearing
  pedagogy (the algebra), the new section is the architectural
  realization.
- Cross-link the partition primitive's consumer table (here:
  `discrete_measures.rst` partition_by) — when the consumer site changes
  semantics ("Wave 2 of SN performance plan; closes Issue #4" →
  "Wave 2 of SN performance plan (closed Issue #4); the SN sweep then
  iterates octants and dispatches each octant to a per-octant
  SweepDependencyGraph — see :ref:`sweep-octant-dependency-graph`"), update
  it. Forward-pointers from primitive sites to consumer sites are the
  glue that lets Nexus see the architecture as a connected graph.
- Trap-fix subsection: name the trap ("L7 trap"), give before/after code
  blocks (legacy `for n in range(N): bc.apply(...)[n]` vs Wave-N
  `for octant in quad.octants: bc.apply(...)`), explain the architectural
  argument WHY the new form is structurally correct (the BC operator's
  semantics are per-octant by construction), AND name the regression-
  detector test with its `@pytest.mark.catches("ERR-NNN")` tag.
- Bit-identity vs principled-equivalence: give a 3-bullet taxonomy when
  applicable (LS-family bit-identical / Lebedev converged-but-not-bit-
  identical / vacuum-BC where partner state doesn't matter). Anchor each
  bullet on a concrete test or verification chain.
- Performance: list-table with ≥3 configs, "Target?" column with honest
  Below/At lower bound/Yes verdicts, "Comment" column explaining the
  regime. Follow with a "follow-up direction" paragraph. NEVER quote a
  speedup without the configuration that produced it.
- Verification chain: L0 unit-test counts (per file), L1 closed-form
  anchor + bug catcher (case names), L2 regression count, with explicit
  reference to the `vv-principles` ladder.
- Don't add `:label:` to subsection HEADERS; reserve `:label:` for
  equations. Use `.. _section-anchor:` for sections (already-existing
  ORPHEUS pattern).
- Add ONE bullet to the Key Facts admonition near the page top so the
  primary architectural change is discoverable from the page header
  (here: "2-D wavefront sweep (Wave 2): per-octant batched dispatch via
  SweepDependencyGraph...").
- Pre-existing Sphinx warnings (the count, not the content): diff before
  vs after the edit. Same count = clean. The line numbers shift due to
  insertions and that's expected; don't try to fix unrelated pre-
  existing warnings during a Wave-N doc update.
- `vv-status: <label> documented` annotation immediately after the
  `.. math:: :label:` block keeps the audit harness happy when the
  equation is a structural framing (eg `streaming-inverse-direct-sum`)
  rather than a tested invariant.
